/**
 * Name: CHEN Qixu
 * Student id: 20412015
 * ITSC email: qchenax@connect.ust.hk
*/
#include <iostream>
#include "mpi_smith_waterman.h"
#include <math.h>

/*
 *  You can add helper functions and variables as you wish.
 */
using namespace std;

int compute_score(int **score, int i, int j, char *a, char *b) {
	return  max(0,
		max(score[i - 1][j - 1] + sub_mat(a[i - 1], b[j - 1]),
			max(score[i - 1][j] - GAP,
				score[i][j - 1] - GAP)));
}

int smith_waterman(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {
	int max_score = 0;
	int local_max_score = 0;

	MPI_Bcast(&a_len, 1, MPI_INT, 0, comm);
	MPI_Bcast(&b_len, 1, MPI_INT, 0, comm);
	MPI_Bcast(&p, 1, MPI_INT, 0, comm);

	if (my_rank != 0) {

		a = new char[a_len + 1];
		b = new char[b_len + 1];
	}
	MPI_Barrier(comm);
	MPI_Bcast(a, a_len + 1, MPI_CHAR, 0, comm);
	MPI_Bcast(b, b_len + 1, MPI_CHAR, 0, comm);

	int total_lv = a_len + b_len;
	int **score = new int*[a_len + 1];
	for (int i = 0; i <= a_len; i++) {
		score[i] = new int[b_len + 1];
		for (int j = 0; j <= b_len; j++) {
			score[i][j] = 0;
		}
	}
	int p_lowerbound = round(my_rank*1.0*(b_len + 1)*1.0 / p);
	int p_upperbound = round((my_rank + 1)*1.0*(b_len + 1)*1.0 / p);
	
	for (int cur_lv = 2; cur_lv <= a_len + b_len; cur_lv++) {
		

		int lowerbound = max(1, cur_lv - a_len);
		int upperbound = min(cur_lv - 1, b_len) + 1;
				
		//cout << "lowerbound for process 0 is " << m_lowerbound << " and upperbound is " << m_upperbound << endl;
		


		int send = 0;
		int recv = 0;
		//need to fetch score[curlv-my_lowerbound][my_lowerbound-1] from previous process
		// correspond to score[curlv-my_upperbound][my_upperbound-1] in the prev process
		if (p_upperbound >= lowerbound && p_lowerbound < upperbound) {
			if (my_rank % 2) {
				if (p_lowerbound >= lowerbound) {
					MPI_Recv(&recv, 1, MPI_INT, my_rank - 1, 0, comm, MPI_STATUS_IGNORE);
					score[cur_lv - p_lowerbound][p_lowerbound - 1] = recv;
				}
				if (p_upperbound < upperbound) {
					send = score[cur_lv - p_upperbound][p_upperbound - 1];
					MPI_Send(&send, 1, MPI_INT, my_rank + 1, 0, comm);
				}
			}
			else {
				if (p_upperbound < upperbound) {
					send = score[cur_lv - p_upperbound][p_upperbound - 1];
					MPI_Send(&send, 1, MPI_INT, my_rank + 1, 0, comm);
				}
				if (p_lowerbound >= lowerbound) {
					MPI_Recv(&recv, 1, MPI_INT, my_rank - 1, 0, comm, MPI_STATUS_IGNORE);
					score[cur_lv - p_lowerbound][p_lowerbound - 1] = recv;
				}
			}
		}


		int my_lowerbound = max(lowerbound, p_lowerbound);
		int my_upperbound = min(upperbound, p_upperbound);
		for (int i = my_lowerbound; i < my_upperbound; i++) {

			score[cur_lv - i][i] = compute_score(score, cur_lv - i, i, a, b);
			local_max_score = max(local_max_score, score[cur_lv - i][i]);
			//cout << "computing for [" << cur_lv - i << "][" << i << "] and score is "<< score[cur_lv - i][i] <<" local max is " << local_max_score<<endl;
		}
		
	}

	MPI_Reduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, 0, comm);
	for (int i = 0; i <= a_len; i++) {
		delete[] score[i];
	}
	delete[] score;

	if (my_rank == 0) {
		return max_score;
	}
	else {
		return 0;
	}
	
}
