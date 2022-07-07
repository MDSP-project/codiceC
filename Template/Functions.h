#include <math.h>
#include "io.h"
#include <fstream>
#include <iostream>
#include "ipp.h"
#include <complex.h>
using namespace std;


#define MAX_FILE_NAME_LENGTH 255


bool read_dat(char* name, double* data, int dim);
void write_dat(char* name, double* data, int dim, char* save_dir);
void cos_h(double* p0, double** H, int M, int N);
void cos_p(double* p0, double** H, int M, int N);
void petr_cos_h(double* p0, double** H, int M, int N);
void analisi(double* InputData_d, double* InputData_x, double** D_buffer, double** X_buffer, double** D, double** X, double** H, double** Hp, int M, int N, int FrameSize);
double* crossfilter(double** X, double** Y, double** X_buffer, double** delay_buffer, double** G, double** G_adj, double** D, int K, int M, int N, int FrameD, int j);
