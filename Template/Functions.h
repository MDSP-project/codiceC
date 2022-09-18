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
void analisi(double* InputData_xl, double* InputData_xr, double** Xl_buffer, double** Xr_buffer, double** Xl, double** Xr, double** H, int M, int N, int FrameSize);
void analisi_petr(double* InputData_x, double** X_buffer, double** X, double** Hp, int M, int N, int FrameSize);
void crossfilter(double** X, double** Y, double** X_buffer, double** W, int K, int M, int N, int FrameD, int j);
void calculatemu(double step_size, Ipp64f* P, double** X,double* mu, int M, double beta, int j);
void adaptation(double** G,double* mu, double* e, double** X_buffer, int K, int M);
void sintesi(double** F, double** Output_Y, double** Y, int M, int N,int Framesize, double* OutputData);
void sintesiE(double** F, double** Output_Y, double** Y, int M, int N, int Framesize, double* OutputData);
void calcG(double** G, double** F, int M, int K, int N);
void filterHRTF(Ipp64f* x, Ipp64f* y, Ipp64f* E_buf, Ipp64f* b, int L, int FrameSize);
void delay(double** X, double** out_delay, double** delay_buffer, double delay_value, int M, int FrameSize);
