// cuda1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#define SEQ_CNT 10
#define SEQ_LEN 10000000 // 10^7

// finds hamming distance between pairs 0-1, 0-2, ..., 0-SEQ_CNT
// and indices of first mismatch
void preprocessing(int *sequences, int *out, int *indices);
// finds hamming disctance between remaining pairs making use of 
// indices and out
void reasoning(int *sequences, int *out, int *indices);
// finds pairs of sequences with hamming distance equals 1
void hamming_distances(int *sequences, int *out);
// process seq1 and seq2 starting from index 1
void hamming(int *sequences, int *out, int *indices, int seq1, int seq2, int start);
// prints bit representation of argument
void bitprint(int x);

int _tmain(int argc, _TCHAR* argv[])
{
	int *sequences = NULL, *out = NULL;

	//====== Memory allocation ======//
	if ((sequences = (int *) malloc(SEQ_CNT * SEQ_LEN * sizeof(int))) == NULL)
	{
		fprintf(stderr, "Out of memory\n");
		getchar();
		exit(EXIT_FAILURE);
	}
	if ((out = (int *) malloc(SEQ_CNT * SEQ_CNT * sizeof(int))) == NULL)
	{
		fprintf(stderr, "Out of memory\n");
		getchar();
		exit(EXIT_FAILURE);
	}
	//===============================//

	//====== Initial values ======//
	srand(17);
	memset(sequences, 0, SEQ_CNT * SEQ_LEN * sizeof(int));
	//for (int k = 0; k < SEQ_CNT; k += 2)
	//sequences[(k + 1) * SEQ_LEN  - 10] = rand()%2 ? 1 : 3;
	for (int k = 0; k < SEQ_CNT; k++)
		sequences[(k + 1) * SEQ_LEN - 100] = rand();
	memset(out, 0, SEQ_CNT * SEQ_CNT * sizeof(int));
	//============================//

	hamming_distances(sequences, out);

	//====== Result ======//
	for (int i = 0; i < SEQ_CNT; i++)
		for (int j = i + 1; j < SEQ_CNT; j++)
			printf("out[%d][%d] = %d\n", i, j, out[i * SEQ_CNT + j]);
	//====================//

	//====== Cleanup ======//
	free(sequences);
	free(out);
	//=====================//

	getchar();
	return 0;
}

void hamming_distances(int *sequences, int *out)
{
	int *indices = NULL;

	if ((indices = (int *) malloc(SEQ_CNT * SEQ_CNT * sizeof(int))) == NULL)
	{
		fprintf(stderr, "Out of memory\n");
		getchar();
		exit(EXIT_FAILURE);
	}

	for (int q = 0; q < SEQ_CNT * SEQ_CNT; q++)
		indices[q] = SEQ_LEN;

	preprocessing(sequences, out, indices);

	reasoning(sequences, out, indices);

	free(indices);
}

void preprocessing(int *sequences, int *out, int *indices)
{
	for (int j = 1; j < SEQ_CNT; j++)
		hamming(sequences, out, indices, 0, j, 0);
}

void reasoning(int *sequences, int *out, int *indices)
{
	for (int i = 1; i < SEQ_CNT; i++)
		for (int j = i + 1; j < SEQ_CNT; j++)
		{
			// i == (i - 1)
			if (out[(i - 1) * SEQ_CNT + i] == 0)
			{
				out[i * SEQ_CNT + j] = out[(i - 1) * SEQ_CNT + j];
				indices[i * SEQ_CNT + j] = indices[(i - 1) * SEQ_CNT + j];
			}
			// j == (i - 1)
			else if (out[(i - 1) * SEQ_CNT + j] == 0)
			{
				out[i * SEQ_CNT + j] = out[(i - 1) * SEQ_CNT + i];
				indices[i * SEQ_CNT + j] = indices[(i - 1) * SEQ_CNT + i];
			}
			// h((i - 1), i) == h((i - 1), j) == 1
			else if (out[(i - 1) * SEQ_CNT + i] == 1 && out[(i - 1) * SEQ_CNT + j] == 1)
			{
				// checking if the mismatch indices are the same
				if (indices[(i - 1) * SEQ_CNT + i] != indices[(i - 1) * SEQ_CNT + j])
				{
					out[i * SEQ_CNT + j] = 2;
					indices[i * SEQ_CNT + j] = std::min(indices[(i - 1) * SEQ_CNT + i], indices[(i - 1) * SEQ_CNT + j]);
				}
				else
				{
					// indices are the same but they correspond to entire integers
					// and we demand every bit to be the same
					int ind = indices[(i - 1) * SEQ_CNT + i];
					unsigned int tmp = sequences[i * SEQ_LEN + ind] ^ sequences[j * SEQ_LEN + ind];
					if (tmp)
					{
						out[i * SEQ_CNT + j] = 2;
						indices[i * SEQ_CNT + j] = ind;
					}
					else
						out[i * SEQ_CNT + j] = 0;
				}
			}
			// (h((i - 1), i) == 1 && h((i - 1), j) >= 2) ||
			// (h((i - 1), i) >= 2 && h((i - 1), j) == 1)
			else if ((out[(i - 1) * SEQ_CNT + i] == 1 && out[(i - 1) * SEQ_CNT + j] >= 2) ||
				(out[(i - 1) * SEQ_CNT + i] >= 2 && out[(i - 1) * SEQ_CNT + j] == 1))
			{
				int ind = std::min(indices[(i - 1) * SEQ_CNT + i], indices[(i - 1) * SEQ_CNT + j]);
				if (indices[(i - 1) * SEQ_CNT + i] != indices[(i - 1) * SEQ_CNT + j])
				{
					out[i * SEQ_CNT + j] = 2;
					indices[i * SEQ_CNT + j] = ind;
				}
				else
					hamming(sequences, out, indices, i, j, ind);
			}
			// (h((i - 1), i) >= 2 && h((i - 1), j) >= 2)
			else
			{
				int ind = std::min(indices[(i - 1) * SEQ_CNT + i], indices[(i - 1) * SEQ_CNT + j]);
				hamming(sequences, out, indices, i, j, ind);
			}
		}
}

void hamming(int *sequences, int *out, int *indices, int seq1, int seq2, int start)
{
	unsigned int tmp = 0;
	unsigned int mask = 1 << ((sizeof(int) * 8) - 1);

	for (int k = start; k < SEQ_LEN; k++)
		if ((tmp = (sequences[seq1 * SEQ_LEN + k] ^ sequences[seq2 * SEQ_LEN + k])))
		{
			// mismatch for sequences i and j found at position k
			indices[seq1 * SEQ_CNT + seq2] = std::min(indices[seq1 * SEQ_CNT + seq2], k);

			// we must check bit by bit how much they mismatch
			while (mask)
			{
				if (tmp & mask)
					if ((++out[seq1 * SEQ_CNT + seq2]) > 1)
						return;

				mask >>= 1;
			}
			// restore mask
			mask = 1 << ((sizeof(int) * 8) - 1);
		}
}

void bitprint(int x)
{
	unsigned int mask = 1 << ((sizeof(int) * 8) - 1);
	while (mask)
	{
		printf("%d", x & mask ? 1 : 0);
		mask >>= 1;
	}
	printf("\n");
}