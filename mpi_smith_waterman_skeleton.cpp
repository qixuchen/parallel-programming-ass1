/**
 * Name:
 * Student id:
 * ITSC email:
*/
#include <iostream>
#include "mpi_smith_waterman.h"

/*
 *  You can add helper functions and variables as you wish.
 */
using namespace std;

int smith_waterman(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {

	int **score = new int*[a_len + 1];
	for (int i = 0; i <= a_len; i++) {
		score[i] = new int[b_len + 1];
		for (int j = 0; j <= b_len; j++) {
			score[i][j] = 0;
		}
	}
	// main loop
	int max_score = 0;
	for (int i = 1; i <= a_len; i++) {
		for (int j = 1; j <= b_len; j++) {
			score[i][j] = max(0,
				max(score[i - 1][j - 1] + sub_mat(a[i - 1], b[j - 1]),
					max(score[i - 1][j] - GAP,
						score[i][j - 1] - GAP)));
			max_score = max(max_score, score[i][j]);
		}

	}
	for (int i = 0; i <= a_len; i++) {
		delete[] score[i];
	}
	delete[] score;

	return max_score;
}
