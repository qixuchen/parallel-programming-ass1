/**
 * Name:
 * Student id:
 * ITSC email:
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

	
	int *sendcount = new int[p];
	int *displs = new int[p];

	
	for (int cur_lv = 2; cur_lv <= a_len + b_len; cur_lv++) {
		

		int lowerbound = max(1, cur_lv - a_len);
		int upperbound = min(cur_lv - 1, b_len) + 1;
		
		int *res_buffer = new int[upperbound - lowerbound];
		
		//cout << "lowerbound for process 0 is " << m_lowerbound << " and upperbound is " << m_upperbound << endl;
		for (int i = 0; i < p; i++) {
			int p_lowerbound = lowerbound + round(i*1.0*(upperbound - lowerbound)*1.0 / p);
			int p_upperbound = lowerbound + round((i + 1)*1.0*(upperbound - lowerbound)*1.0 / p);
			sendcount[i] = p_upperbound - p_lowerbound;
			displs[i] = p_lowerbound - lowerbound;
			//cout << "lowerbound for process " << i << " is " << p_lowerbound << " and upperbound is " << p_upperbound << endl;
		}
		
		int my_lowerbound = lowerbound + round(my_rank*1.0*(upperbound - lowerbound)*1.0 / p);
		int my_upperbound = lowerbound + round((my_rank + 1)*1.0*(upperbound - lowerbound)*1.0 / p);
		int *my_part = new int[my_upperbound - my_lowerbound];
		for (int i = my_lowerbound; i < my_upperbound; i++) {
			
			score[cur_lv - i][i] = compute_score(score, cur_lv - i, i, a, b);
			my_part[i - my_lowerbound] = score[cur_lv - i][i];
			local_max_score = max(local_max_score, score[cur_lv - i][i]);
			//cout << "computing for [" << cur_lv - i << "][" << i << "] and score is "<< score[cur_lv - i][i] <<" local max is " << local_max_score<<endl;
		}
		//cout << "localmax for process " << my_rank << " is " << local_max_score << " and max is " << max_score << endl;*/

		
		
		MPI_Barrier(comm);
		MPI_Allgatherv(my_part, sendcount[my_rank], MPI_INT, res_buffer, sendcount, displs, MPI_INT, comm);
		
		for (int i = 0; i < upperbound - lowerbound; i++) {
			score[cur_lv - lowerbound - i][lowerbound + i] = res_buffer[i];
		}
		

		delete[] my_part;
		my_part = nullptr;
		delete[] res_buffer;
		res_buffer = nullptr;
		
		

	}

	MPI_Reduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, 0, comm);
	delete[] sendcount;
	sendcount = nullptr;
	delete[] displs;
	displs = nullptr;
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
