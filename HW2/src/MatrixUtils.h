#pragma once

int index(int i, int j, int m);
void generate_random_matrix(double* M, int n, int m);
void generate_random_convergent_matrix(double* M, int n, int m);
void print_matrix(double* M, int n, int m);
void swap_rows(double* M, int m, int i, int k);
double rand_coeff();
double* transpose(double* M, int n, int m);