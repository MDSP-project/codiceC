#include "stdafx.h"
#include "Functions.h"


void write_dat(char* name, double* data, int dim, char* save_name) {

	char name_a[MAX_FILE_NAME_LENGTH];
	memset(name_a, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	strcpy(name_a, save_name);
	strcat(name_a, name);

	char* c;
	c = (char*)(void*)data;

	fstream File;
	File.open(name_a, ios::out | ios::binary);
	File.write(c, dim * sizeof(double));
	File.close();
}

bool read_dat(char* name, double* data, int dim) {
	char* c;
	c = (char*)(void*)data;

	fstream File;
	bool flag = false;
	File.open(name, ios::in | ios::binary);
	if (File.is_open())
	{
		flag = true;
	}
	File.read(c, dim * sizeof(double));
	File.close();
	return flag;
}


void cos_h(double* p0, double** H, int M, int N) {

	double* tmp;
	tmp = ippsMalloc_64f(N);

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {

			tmp[j] = 2.0 * cos((IPP_PI / ((double)M)) * ((double)i + 0.5) * ((double)j - (((double)N - 1.0) / 2.0)) + (pow(-1.0, (double)i)) * (IPP_PI / 4.0));
		}
		ippsMul_64f(p0, tmp, H[i], N);
	}


	ippsFree(tmp);

}

void cos_p(double* p0, double** H, int M, int N) {

	double* tmp;
	tmp = ippsMalloc_64f(N);


	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {

			tmp[j] = 2.0 * cos((IPP_PI / ((double)M)) * ((double)i + 0.5) * ((double)j - (((double)N - 1.0) / 2.0)) - (pow(-1.0, (double)i)) * (IPP_PI / 4.0));
		}
		ippsMul_64f(p0, tmp, H[i], N);
	}


	ippsFree(tmp);

}


void petr_cos_h(double* p0, double** H, int M, int N) {

	/*
	Questa funzione effettua la modulazione coseno del filtro prototipo p0,
	per generare il banco di analisi secondo lo schema della Petraglia.

	M � il numero di bande
	N � la dimensione nel tempo
	*/
	Ipp64f* q1_re, *q2_re, *q1_im, *q2_im, * q0_tmp, *p0_tmp, *q00, *q01, *phase1, *phase2;
	Ipp8u* pBuffer;
	int bufSize = 0;
	IppStatus status;
	IppEnum funCfg = (IppEnum)(ippAlgAuto);
	double* tmp_1, *tmp_2;
	

	tmp_1 = ippsMalloc_64f(2 * N - 1);
	tmp_2 = ippsMalloc_64f(2 * N - 1);
	q1_re = ippsMalloc_64f(2 * N);
	q1_im = ippsMalloc_64f(2 * N);
	q2_re = ippsMalloc_64f(2 * N);
	q2_im = ippsMalloc_64f(2 * N);
	q0_tmp = ippsMalloc_64f(2 * N);
	p0_tmp = ippsMalloc_64f(2 * N);
	q00 = ippsMalloc_64f(2 * N);
	q01 = ippsMalloc_64f(2 * N);
	phase1 = ippsMalloc_64f(N);
	phase2 = ippsMalloc_64f(N);


	ippsZero_64f(q1_re , 2 * N);
	ippsZero_64f(q1_im, 2 * N);
	ippsZero_64f(q2_re, 2 * N);
	ippsZero_64f(q2_im, 2 * N);
	ippsZero_64f(q0_tmp, 2 * N);



	for (int n = 0; n < N; n++) {
		phase1[n] = - (IPP_PI / (2 * (double)M))      * (double)n ;
		phase2[n] = (IPP_PI / (2 * (double)M)) * (double)n;
	}

	ippsPolarToCart_64f(p0, phase1, q1_re, q1_im, N);
	ippsPolarToCart_64f(p0, phase2, q2_re, q2_im, N);

	// convoluzione troncata ai reali (percorsi incrociati)
	for (int i = 0; i < 2*N - 1; i++) {
		for (int j = 0; j < N; j++) {
			if ((i >= j) && (i - j) < N) {

				q0_tmp[i] = q0_tmp[i] + q1_re[j] * q2_re[i - j] - q1_im[j] * q2_im[i - j]; 

			}

		}
	}

	status = ippsConvolveGetBufferSize(N,N,ipp64f,funCfg,&bufSize);  

	pBuffer = ippsMalloc_8u(bufSize);

	ippsConvolve_64f(p0,N,p0,N,p0_tmp,0,pBuffer);  //convoluzione per i percorsi dritti (numeri reali)

	int q = 0;

	for (int m = 0; m < M; m++) {
		
		for (int n = 0; n < (2*N -1); n++) {
			tmp_1[n] = 2.0 * cos((IPP_PI / ((double)M)) * ((double)m + 0.5) * ((double)n - ((double)N / 2.0)) + 2*((pow(-1.0, (double)m)) * (IPP_PI / 4.0)));
			
			if (m < M-1) tmp_2[n] = 2.0 * cos((IPP_PI / ((double)M)) * ((double)m + 1) * (((double)n) - ((double)N / 2)));

			
		}

		
		ippsMul_64f(p0_tmp, tmp_1, H[q], 2*N - 1);
		if (q < (2*M-2)) ippsMul_64f(q0_tmp, tmp_2, H[q + 1], 2*N-1);
		q = q + 2;

		
	}

	/*
	Per k = 0 e k = M-1 bisogna considerare dei termini aggiuntivi come
    indicato nel paper:
    "New Structures for Adaptive Filtering in Subbands with Critical Sampling" 
    by Mariane R. Petraglia
	*/

	for (int n = 0; n < 2 * M -1; n++) {


		q00[n] = 2 * q0_tmp[n];
		q01[n] = ((double)M - 2) * (double)n * q0_tmp[n];
		H[0][n] = H[0][n] + q00[n];
		H[2*M-2][n] = H[2 * M - 2][n] + q01[n];




	}

	// sezione di disallocazione vettore temporanei
	
	ippsFree(tmp_1);
    ippsFree(tmp_2);
	ippsFree(q1_re);
	ippsFree(q1_im);
	ippsFree(q2_re);
	ippsFree(q2_im);
	ippsFree(q0_tmp);
	ippsFree(p0_tmp);
	ippsFree(q00);
	ippsFree(q01);
 	ippsFree(phase1);
	ippsFree(phase2);

}



void analisi(double* InputData_xl, double* InputData_xr, double** Xl_buffer, double** Xr_buffer, double** Xl, double** Xr, double** H, int M, int N, int FrameSize) {
	
	Ipp64f* xl = 0;
	Ipp64f* xr = 0;
	double* temp1 = 0;
	double* temp2 = 0;
	double* U1 = 0;
	double* U2 = 0;
	int i = 0;



	if (xl == 0)
	{
		xl = ippsMalloc_64f(FrameSize);
		ippsZero_64f(xl, FrameSize);
	}

	if (xr == 0)
	{
		xr = ippsMalloc_64f(FrameSize);
		ippsZero_64f(xr, FrameSize);
	}

	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(N);
		ippsZero_64f(temp1, N);
	}

	if (temp2 == 0)
	{
		temp2 = ippsMalloc_64f(N);
		ippsZero_64f(temp2, N);
	}

	if (U1 == 0)
	{
		U1 = ippsMalloc_64f(M);
		ippsZero_64f(U1, M);
	}

	if (U2 == 0)
	{
		U2 = ippsMalloc_64f(M);
		ippsZero_64f(U2, M);
	}

	for (int n = 0; n < FrameSize; n++)
	{
		xl[n] = (double)InputData_xl[n] / 32768.0;
		xr[n] = (double)InputData_xr[n] / 32768.0;

		for (int m = 0; m < M ; m++)
		{
		
			// scorrimento del vettore X_buff e aggiornamento
			ippsCopy_64f(Xl_buffer[m], temp1, N);
			ippsCopy_64f(temp1, Xl_buffer[m] + 1, N);
			Xl_buffer[m][0] = xl[n];

			// moltiplicazione per il banco di analisi
			//ippsDotProd_64f(H[m], X_buff[m], N, &Hx[m]);
			ippsMul_64f(H[m], Xl_buffer[m], temp1, N);
			ippsSum_64f(temp1, N, &U1[m]);

			// scorrimento del vettore X_buff e aggiornamento
			ippsCopy_64f(Xr_buffer[m], temp2, N);
			ippsCopy_64f(temp1, Xr_buffer[m] + 1, N);
			Xr_buffer[m][0] = xr[n];

			// moltiplicazione per il banco di analisi
			//ippsDotProd_64f(H[m], X_buff[m], N, &Hx[m]);
			ippsMul_64f(H[m], Xr_buffer[m], temp2, N);
			ippsSum_64f(temp2, N, &U2[m]);
		

			//decimazione e interpolazione
		if ((n % M) == 0) 
			{
			Xl[m][i] = U1[m]; //downsampling
			Xr[m][i] = U2[m]; //downsampling
		    }

		
		}
		
		if ((n % M) == 0)
		{
			i++; //downsampling

		}
	}


	if (xl != 0)
	{
		ippsFree(xl);
		xl = 0;
	}

	if (xr != 0)
	{
		ippsFree(xr);
		xr = 0;
	}

	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}

	if (temp2 != 0)
	{
		ippsFree(temp2);
		temp2 = 0;
	}

	if (U1 != 0)
	{
		ippsFree(U1);
		U1 = 0;
	}
	
	if (U2 != 0)
	{
		ippsFree(U2);
		U2 = 0;
	}

}

void analisi_petr(double* InputData_x, double** X_buffer, double** X, double** Hp, int M, int N, int FrameSize) {

	Ipp64f* x = 0;
	double* temp1 = 0;
	double* U1 = 0;
	int i = 0;



	if (x == 0)
	{
		x = ippsMalloc_64f(FrameSize);
		ippsZero_64f(x, FrameSize);
	}

	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(2*N-1);
		ippsZero_64f(temp1, 2*N-1);
	}

	if (U1 == 0)
	{
		U1 = ippsMalloc_64f(2*M-1);
		ippsZero_64f(U1, 2*M-1);
	}

	for (int n = 0; n < FrameSize; n++)
	{
		x[n] = (double)InputData_x[n] / 32768.0;
	

		for (int m = 0; m < 2*M-1; m++)
		{

			// scorrimento del vettore X_buff e aggiornamento
			ippsCopy_64f(X_buffer[m], temp1, 2*N-1);
			ippsCopy_64f(temp1, X_buffer[m] + 1, 2*N-2);
			X_buffer[m][0] = x[n];

			// moltiplicazione per il banco di analisi
			//ippsDotProd_64f(H[m], X_buff[m], N, &Hx[m]);
			ippsMul_64f(Hp[m], X_buffer[m], temp1, 2*N-1);
			ippsSum_64f(temp1, 2*N-1, &U1[m]);



			//decimazione e interpolazione
			if ((n % M) == 0)
			{
				X[m][i] = U1[m]; //downsampling
			}


		}

		if ((n % M) == 0)
		{
			i++; //downsampling

		}
	}


	if (x != 0)
	{
		ippsFree(x);
		x = 0;
	}

	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}

	if (U1 != 0)
	{
		ippsFree(U1);
		U1 = 0;
	}

}

void crossfilter(double** X, double** Y, double** X_buffer, double** W, int K, int M, int N, int FrameD, int j)
{
	Ipp64f* temp1=0;
	Ipp64f* Y_tmp = 0;

	

	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(K);
		ippsZero_64f(temp1, K);
	}

	if (Y_tmp == 0)
	{
		Y_tmp = ippsMalloc_64f(M);
		ippsZero_64f(Y_tmp, M);
	}


	for (int m = 0; m < 2 * M - 1; m++)
	{
		// scorrimento del vettore X_buff e aggiornamento
		ippsCopy_64f(X_buffer[m], temp1,K);
		ippsCopy_64f(temp1, X_buffer[m] + 1, K - 1); 
		X_buffer[m][0] = X[m][j];
	}

	ippsZero_64f(Y_tmp, M);

	int q = 2;

	for (int m = 1; m < M - 1; m++)
	{
		for(int k=0; k<K; k++)
		{
			Y_tmp[m] = Y_tmp[m] + X_buffer[q][k] * W[m][k] + X_buffer[q - 1][k] * W[m - 1][k] + X_buffer[q + 1][k] * W[q + 1][k];
		}
		q = q + 2;
	}
	for (int k = 0; k < K; k++)
	{
		Y_tmp[0] = Y_tmp[0] + X_buffer[0][k] * W[0][k] + X_buffer[1][k] * W[1][k];
		Y_tmp[M-1]=Y_tmp[M-1]+ X_buffer[2*M-2][k] * W[M-1][k] + X_buffer[2*M-3][k] * W[M-2][k];
	}
	for (int m = 0; m < M; m++)
	{
		Y[m][j] = Y_tmp[m];
		
	}

	
	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}

	if (Y_tmp != 0)
	{
		ippsFree(Y_tmp);
		Y_tmp = 0;
	}
	
}


void calculatemu(double step_size, Ipp64f* P1, Ipp64f* P2, double** X1, double** X2, Ipp64f* mu, int M, double alpha, double beta, int j)
{
	for (int m = 0; m < 2*M-1; m++)
	{
		P1[m] = P1[m] * beta + (1 - beta) * (abs(pow(X1[m][j], 2)));
		P2[m] = P2[m] * beta + (1 - beta) * (abs(pow(X2[m][j], 2)));
	}

	mu[0] = step_size / (alpha + P1[0] + P2[0] + P1[1] + P2[1]);
	mu[2 * M - 2] = step_size / (alpha + P1[2 * M - 2] + P1[2 * M - 3] + P2[2 * M - 2] + P2[2 * M - 3]);
	for (int m = 1; m < 2*M-2; m++)
	{
		mu[m] = step_size / (alpha + P1[m - 1] + P1[m] + P1[m + 1] + P2[m - 1] + P2[m] + P2[m + 1]);
	}
}

void adaptation(double** W, double* mu, double* e1, double* e2, double** X1, double** X2, int K, int M)
{
	double** W_adj = 0;
	int q = 2;
	W_adj = new double* [M];
	for (int i = 0; i < M; i++)
	{
		W_adj[i] = new double[K];
		memset(W_adj[i], 0.0, (K) * sizeof(double));
	}

	for (int m = 1; m < M-1; m++)
	{
		for (int k = 0; k < K; k++)
		{
			W_adj[m][k] = W[m][k] + mu[m] * (e1[m] * X1[q][k] + e1[m - 1] * X1[q - 1][k] +  e1[m + 1] * X1[q + 1][k] + e2[m] * X2[q][k] + e2[m - 1] * X2[q - 1][k] + e2[m + 1] * X2[q + 1][k]);
		}
		q = q + 2;
	}

	for (int k = 0; k < K; k++)
	{
		W_adj[0][k] = W[0][k] + mu[0] * (e1[0] * X1[0][k] + e1[1] * X1[1][k]+ e2[0] * X2[0][k] + e2[1] * X2[1][k]);
		W_adj[M-1][k] = W[M-1][k] + mu[M-1] * (e1[M-1] * X1[M-1][k] + e1[M-2] * X1[M-2][k]+ e2[M - 1] * X2[M - 1][k] + e2[M - 2] * X2[M - 2][k]);
	}

	for (int m = 0; m < M; m++)
	{
		ippsCopy_64f(W_adj[m], W[m], K);
	}

	for (int i = 0; i < M; i++)
		delete[] W_adj[i];
	delete[] W_adj;
}

void sintesi(double** F, double** Output_Y, double** Y, int M, int N,int Framesize, double* OutputData)
{
	double** interp = 0;
	int i = 0;
	Ipp64f* temp1= 0;
	Ipp64f* Gw = 0;
	Ipp64f* y = 0;

	interp = new double* [M];
	for (int i = 0; i < M; i++)
	{
		interp[i] = new double[Framesize];
		memset(interp[i], 0.0, (Framesize) * sizeof(double));
	}


	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(N);
		ippsZero_64f(temp1, N);
	}

	if (Gw == 0)
	{
		Gw = ippsMalloc_64f(M);
		ippsZero_64f(Gw, M);
	}

	if (y == 0)
	{
		y = ippsMalloc_64f(Framesize);
		ippsZero_64f(y, Framesize);
	}

	for (int n = 0; n < Framesize; n++)
	{
		for (int m = 0; m < M; m++)
		{
			if ((n % M) == 0) 
			{
				interp[m][n] = Y[m][i]; //downsampling
				if (m==15)
				{
					i = i + 1;
				}
			}
			else
				interp[m][n] = 0; //upsampling

			
		}

		for (int m = 0; m < M; m++)
		{
			ippsCopy_64f(Output_Y[m], temp1, N);
			ippsCopy_64f(temp1, Output_Y[m] + 1, N - 1);
			Output_Y[m][0] = interp[m][n];

			ippsDotProd_64f(F[m], Output_Y[m], N, &Gw[m]);
		}

		ippsSum_64f(Gw, M, &y[n]);
		OutputData[n] = y[n] * 32768.0 * M;
		
	}	

	for (int i = 0; i < M; i++)
		delete[] interp[i];
	delete[] interp;

	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}

	if (y != 0)
	{
		ippsFree(y);
		y = 0;
	}

	if (Gw != 0)
	{
		ippsFree(Gw);
		Gw = 0;
	}
}

void sintesiE(double** F, double** Output_Y, double** Y, int M, int N, int Framesize, double* OutputData)
{
	double** interp = 0;
	int i = 0;
	Ipp64f* temp1 = 0;
	Ipp64f* Gw = 0;
	Ipp64f* y = 0;

	interp = new double* [M];
	for (int i = 0; i < M; i++)
	{
		interp[i] = new double[Framesize];
		memset(interp[i], 0.0, (Framesize) * sizeof(double));
	}


	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(N);
		ippsZero_64f(temp1, N);
	}

	if (Gw == 0)
	{
		Gw = ippsMalloc_64f(M);
		ippsZero_64f(Gw, M);
	}

	if (y == 0)
	{
		y = ippsMalloc_64f(Framesize);
		ippsZero_64f(y, Framesize);
	}

	for (int n = 0; n < Framesize; n++)
	{
		for (int m = 0; m < M; m++)
		{
			if ((n % M) == 0)
			{
				interp[m][n] = Y[m][i]; //downsampling
				if (m == 15)
				{
					i = i + 1;
				}
			}
			else
				interp[m][n] = 0; //upsampling


		}

		for (int m = 0; m < M; m++)
		{
			ippsCopy_64f(Output_Y[m], temp1, N);
			ippsCopy_64f(temp1, Output_Y[m] + 1, N - 1);
			Output_Y[m][0] = interp[m][n];

			ippsDotProd_64f(F[m], Output_Y[m], N, &Gw[m]);
		}

		ippsSum_64f(Gw, M, &y[n]);
		OutputData[n] = y[n];

	}

	for (int i = 0; i < M; i++)
		delete[] interp[i];
	delete[] interp;

	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}

	if (y != 0)
	{
		ippsFree(y);
		y = 0;
	}

	if (Gw != 0)
	{
		ippsFree(Gw);
		Gw = 0;
	}
}


void calcG( double** G, double** F, int M, int K, int N)
{

	Ipp8u* pBuffer;
	int bufSize = 0;
	IppStatus status;
	IppEnum funCfg = (IppEnum)(ippAlgAuto);

	Ipp64f* g;
	g = ippsMalloc_64f(K + N - 1);
	ippsZero_64f(g, K + N - 1);

	Ipp64f* g_tmp;
	g_tmp = ippsMalloc_64f(K+N-1);
	ippsZero_64f(g_tmp, K+N-1);

	for (int m = 0; m < M; m++)
	{
		status = ippsConvolveGetBufferSize(K, N, ipp64f, funCfg, &bufSize);
		pBuffer = ippsMalloc_8u(bufSize);
		ippsConvolve_64f(G[m], K, F[m], N, g_tmp, 0, pBuffer);  //convoluzione per i percorsi dritti (numeri reali)

		for (int k = 0; k < K+N-1; k++)
		{
			g[k] += g_tmp[k];

		}

	}

	write_dat("\\riposta.dat", g, K + N - 1, "C:\\Users\\alleg\\Desktop\\Nuova cartella");

	ippsFree(g_tmp);
	ippsFree(g);
}

void filterHRTF(Ipp64f* x, Ipp64f* y, Ipp64f* E_buf, Ipp64f* b, int L, int FrameSize)
{
	double* temp1 = 0;

	if (temp1 == 0)
	{
		temp1 = ippsMalloc_64f(L);
		ippsZero_64f(temp1, L);
	}

	for (size_t i = 0; i < FrameSize; i++)
	{
		double d = 0.0;

		ippsCopy_64f(E_buf, temp1, L);
		ippsCopy_64f(temp1, E_buf + 1, L);
		E_buf[0] = x[i];

		for (size_t n = 0; n < L; n++)
		{

			d = d+E_buf[n]*b[n];

		}

		y[i] = d;

	}

	if (temp1 != 0)
	{
		ippsFree(temp1);
		temp1 = 0;
	}
}

void delay(double** X, double** out_delay, double** delay_buffer, double delay_value, int M, int FrameSize)
{
	for (size_t m = 0; m < M; m++)
	{
		if (delay_value < FrameSize)
		{

			ippsCopy_64f(delay_buffer[m], out_delay[m], delay_value);
			
		}
		else
		{
			ippsCopy_64f(delay_buffer[m], out_delay[m], FrameSize);
		}

		if (delay_value > FrameSize)
		{

			ippsCopy_64f(delay_buffer[m], out_delay[m]+FrameSize+1, delay_value-FrameSize);
			ippsCopy_64f(delay_buffer[m] + int(delay_value) - FrameSize + 1, X[m], FrameSize);
		}
		else if(delay_value!=FrameSize )
		{
			ippsCopy_64f(out_delay[m] + int(delay_value), X[m], FrameSize - delay_value);
		}
		ippsCopy_64f(delay_buffer[m], X[m] + FrameSize - int(delay_value), FrameSize);

	}
}

