#include "StdAfx.h"
#include ".\plugin.h"


//FINALE

PlugIn::PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg): LEEffect(_CBFunction,_PlugRef,ParentDlg)
{

	LESetNumInput(2);  //dichiarazione 2 ingressi
	LESetNumOutput(4); //dichiarazione 2 uscite

	FrameSize = CBFunction(this,NUTS_GET_FS_SR,0,(LPVOID)AUDIOPROC);
	SampleRate = CBFunction(this,NUTS_GET_FS_SR,1,(LPVOID)AUDIOPROC);	

	N = 256;  // lunghezza filtro prototipo
	M = 16;   // numero Bande
	L = 1024;  // lunghezza filtro incognito
	step_size = 0.0005; //Valore massimo dello step size per l'adattamento
	beta = 0.99; // Peso per la stima della potenza in ogni banda
	alpha = 0.5E-10;
	Buf_dim = FrameSize;
	K =256; // Numero di tappi per ogni filtro adattivo
	delay_value = N / M; // Valore del ritardo

	p0 = 0;
	i = 0;
	hll, hlr, hrl, hrr, hrtf, hrtfmax = 0;
	h11xl, h21xl, h12xl, h22xl = 0;
	h11xr, h21xr, h12xr, h22xr = 0;
	h11_sx_buffer, h21_sx_buffer, h12_sx_buffer, h22_sx_buffer = 0;
	h11_ds_buffer, h21_ds_buffer, h12_ds_buffer, h22_ds_buffer = 0;
	P1_1, P2_1, P3_1, P4_1 = 0;
	P1_2, P2_2, P3_2, P4_2 = 0;
	mu1, mu2, mu3, mu4 = 0;
	e1, e2 = 0;
	y1, y2 = 0;
	error1, error2 = 0;

	memset(save_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	memset(hll_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	memset(hlr_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	memset(hrl_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	memset(hrr_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));

	strcpy(save_name, "C:\\Users\\alleg\\Desktop\\SACC_Codice Matlab 14 Luglio\\prototipoMC.dat");
	strcpy(hll_name, "C:\\Users\\alleg\Desktop\\SACC\\p1_left_ch0.f64");
	strcpy(hlr_name, "C:\\Users\\alleg\Desktop\\SACC\\p1_left_ch1.f64");
	strcpy(hrl_name, "C:\\Users\\alleg\Desktop\\SACC\\p1_right_ch0.f64");
	strcpy(hrr_name, "C:\\Users\\alleg\Desktop\\SACC\\p1_right_ch1.f64");
}

int __stdcall PlugIn::LEPlugin_Process(PinType **Input,PinType **Output,LPVOID ExtraInfo)
{ 
	
	double* InputData_xl= ((double*)Input[0]->DataBuffer);
	double* InputData_xr = ((double*)Input[1]->DataBuffer);

	double* OutputDataL = ((double*)Output[0]->DataBuffer);
	double* OutputDataR = ((double*)Output[1]->DataBuffer); 
	double* OutputDataE_L = ((double*)Output[2]->DataBuffer);
	double* OutputDataE_R = ((double*)Output[3]->DataBuffer);


	analisi(InputData_xl, InputData_xr, an_buffer1, an_buffer2, X1, X2, H, M,  N,  FrameSize);   //Funzione blocchi di analisi

	//------------------------------------------------------------- HRTF
	filterHRTF(InputData_xl, h11xl, h11_sx_buffer, hll, L,FrameSize);
	filterHRTF(InputData_xl, h12xl, h12_sx_buffer, hlr, L, FrameSize);
	filterHRTF(InputData_xl, h21xl, h21_sx_buffer, hrl, L, FrameSize);
	filterHRTF(InputData_xl, h22xl, h22_sx_buffer, hrr, L, FrameSize);

	filterHRTF(InputData_xr, h11xr, h11_ds_buffer, hll, L, FrameSize);
	filterHRTF(InputData_xr, h12xr, h12_ds_buffer, hlr, L, FrameSize);
	filterHRTF(InputData_xr, h21xr, h21_ds_buffer, hrl, L, FrameSize);
	filterHRTF(InputData_xr, h22xr, h22_ds_buffer, hrr, L, FrameSize);

	//----------------------------------------------------------------  PETRAGLIA
	analisi_petr(h11xl, h11xl_buffer, out_M1_1, Hp, M, N, FrameSize);
	analisi_petr(h12xl, h12xl_buffer, out_M2_1, Hp, M, N, FrameSize);
	analisi_petr(h11xr, h11xr_buffer, out_M3_1, Hp, M, N, FrameSize);
	analisi_petr(h12xr, h12xr_buffer, out_M4_1, Hp, M, N, FrameSize);

	analisi_petr(h22xr, h22xr_buffer, out_M1_2, Hp, M, N, FrameSize);
	analisi_petr(h21xr, h21xr_buffer, out_M2_2, Hp, M, N, FrameSize);
	analisi_petr(h22xl, h22xl_buffer, out_M3_2, Hp, M, N, FrameSize);
	analisi_petr(h21xl, h21xl_buffer, out_M4_2, Hp, M, N, FrameSize);

	//---------------------------------------------------------------- DELAY

	delay(X1, out_buf_dly1, dly1, delay_value, M, FrameSize);
	delay(X2, out_buf_dly2, dly2, delay_value, M, FrameSize);

	//----------------------------------------------------------------

	for (int j = 0; j < FrameD; j++)
	{
		crossfilter(out_M1_1, out_w1_1, Z_w1_1, W1, K, N, N, FrameD, j);
		crossfilter(out_M2_1, out_w2_1, Z_w2_1, W2, K, N, N, FrameD, j);
		crossfilter(out_M3_1, out_w3_1, Z_w3_1, W3, K, N, N, FrameD, j);
		crossfilter(out_M4_1, out_w4_1, Z_w4_1, W4, K, N, N, FrameD, j);

		crossfilter(out_M1_2, out_w4_2, Z_w4_2, W1, K, N, N, FrameD, j);
		crossfilter(out_M2_2, out_w3_2, Z_w3_2, W2, K, N, N, FrameD, j);
		crossfilter(out_M3_2, out_w2_2, Z_w2_2, W3, K, N, N, FrameD, j);
		crossfilter(out_M4_2, out_w1_2, Z_w1_2, W4, K, N, N, FrameD, j);

		if (i>0)
		{
			calculatemu(step_size, P1_1, P1_2, Z_w1_1, Z_w1_2, mu1, M, alpha, beta, j);
			calculatemu(step_size, P2_1, P2_2, Z_w2_1, Z_w2_2, mu2, M, alpha, beta, j);
			calculatemu(step_size, P3_1, P3_2, Z_w3_1, Z_w3_2, mu3, M, alpha, beta, j);
			calculatemu(step_size, P4_1, P4_2, Z_w4_1, Z_w4_2, mu4, M, alpha, beta, j);

			for (size_t m = 0; m < M; m++)
			{
				out_sum_1[m][j] = out_w1_1[m][j] + out_w2_1[m][j] + out_w3_1[m][j] + out_w4_1[m][j];
				out_sum_2[m][j] = out_w1_2[m][j] + out_w2_2[m][j] + out_w3_2[m][j] + out_w4_2[m][j];
				e1[m] = out_buf_dly1[m][j] - out_sum_1[m][j];
				e2[m] = out_buf_dly2[m][j] - out_sum_2[m][j];
				E1[m][j] = e1[m];
				E2[m][j] = e2[m];
			}

			adaptation(W1, mu1, e1, e2, Z_w1_1, Z_w1_2, K, M);
			adaptation(W2, mu2, e1, e2, Z_w2_1, Z_w2_2, K, M);
			adaptation(W3, mu3, e1, e2, Z_w3_1, Z_w3_2, K, M);
			adaptation(W4, mu4, e1, e2, Z_w4_1, Z_w4_2, K, M);

			
		}
		
	} 
	
	sintesi(F, y1_buf, out_sum_1, M, N, FrameSize, OutputDataL);    //Funzione blocchi di sintesi
	sintesi(F, y2_buf, out_sum_2, M, N, FrameSize, OutputDataR);
	sintesiE(F, error1_buf, E1, M, N, FrameSize, OutputDataE_L);
	sintesiE(F, error2_buf, E2, M, N, FrameSize, OutputDataE_R);

	i=1;  // dopo il primo frame adatto i filtri G

	return COMPLETED;
}

void __stdcall PlugIn::LEPlugin_Init()
{
	//Inizializzazione dei Filtri Adattivi e del vettore delle Potenze

	FrameD = FrameSize / M; // Dimensione Frame decimato

	if (p0 == 0) {
		p0 = ippsMalloc_64f(N);
		ippsZero_64f(p0,N);
	}

	 

// inizializzo HRTF

	if (hll == 0) {
		hll = ippsMalloc_64f(L);
		ippsZero_64f(hll, L);
	}

	if (hlr == 0) {
		hlr = ippsMalloc_64f(L);
		ippsZero_64f(hlr, L);
	}

	if (hrl == 0) {
		hrl = ippsMalloc_64f(L);
		ippsZero_64f(hrl, L);
	}

	if (hrr == 0) {
		hrr = ippsMalloc_64f(L);
		ippsZero_64f(hrr, L);
	}

	if (hrtf == 0) {
		hrtf = ippsMalloc_64f(4*L);
		ippsZero_64f(hrtf, 4*L);
	}

	if (hrtfmax == 0) {
		hrtfmax = ippsMalloc_64f(1);
		ippsZero_64f(hrtfmax ,1);
	}

//----------------------------------------

	if (h11xl == 0) {
		h11xl = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h11xl, Buf_dim);
	}

	if (h12xl == 0) {
		h12xl = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h12xl, Buf_dim);
	}

	if (h21xl == 0) {
		h21xl = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h21xl, Buf_dim);
	}

	if (h22xl == 0) {
		h22xl = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h22xl, Buf_dim);
	}

	//--------------------------------------------


	if (h11xr == 0) {
		h11xr = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h11xr, Buf_dim);
	}

	if (h12xr == 0) {
		h12xr = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h12xr, Buf_dim);
	}

	if (h21xr == 0) {
		h21xr = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h21xr, Buf_dim);
	}

	if (h22xr == 0) {
		h22xr = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(h22xr, Buf_dim);
	}

	//--------------------------------------------

	if (h11_ds_buffer == 0) {
		h11_ds_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h11_ds_buffer, L);
	}

	if (h12_ds_buffer == 0) {
		h12_ds_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h12_ds_buffer, L);
	}

	if (h21_ds_buffer == 0) {
		h21_ds_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h21_ds_buffer, L);
	}

	if (h22_ds_buffer == 0) {
		h22_ds_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h22_ds_buffer, L);
	}

	//---------------------------------------------

	if (h11_sx_buffer == 0) {
		h11_sx_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h11_sx_buffer, L);
	}

	if (h12_sx_buffer == 0) {
		h12_sx_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h12_sx_buffer, L);
	}

	if (h21_sx_buffer == 0) {
		h21_sx_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h21_sx_buffer, L);
	}

	if (h22_sx_buffer == 0) {
		h22_sx_buffer = ippsMalloc_64f(L);
		ippsZero_64f(h22_sx_buffer, L);
	}

	//---------------------------------------------

	if (P1_1 == 0) {
		P1_1 = ippsMalloc_64f(2*M-1);
		ippsZero_64f(P1_1, 2*M-1);
	}

	if (P2_1 == 0) {
		P2_1 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P2_1, 2 * M - 1);
	}

	if (P3_1 == 0) {
		P3_1 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P3_1, 2 * M - 1);
	}

	if (P4_1 == 0) {
		P4_1 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P4_1, 2 * M - 1);
	}
	//---------------------------------------------


	if (P1_2 == 0) {
		P1_2 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P1_2, 2 * M - 1);
	}

	if (P2_2 == 0) {
		P2_2 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P2_2, 2 * M - 1);
	}

	if (P3_2 == 0) {
		P3_2 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P3_2, 2 * M - 1);
	}

	if (P4_2 == 0) {
		P4_2 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(P4_2, 2 * M - 1);
	}
	//---------------------------------------------

	if (mu1 == 0) {
		mu1 = ippsMalloc_64f(M);
		ippsZero_64f(mu1, M);
	}

	if (mu2 == 0) {
		mu2 = ippsMalloc_64f(M);
		ippsZero_64f(mu2, M);
	}

	if (mu3 == 0) {
		mu3 = ippsMalloc_64f(M);
		ippsZero_64f(mu3, M);
	}

	if (mu4 == 0) {
		mu4 = ippsMalloc_64f(M);
		ippsZero_64f(mu4, M);
	}

	//---------------------------------------------

	if (e1 == 0) {
		e1 = ippsMalloc_64f(M);
		ippsZero_64f(e1, M);
	}

	if (e2 == 0) {
		e2 = ippsMalloc_64f(M);
		ippsZero_64f(e2, M);
	}

	//------------------------------------------------

	if (y1 == 0) {
		y1 = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(y1, Buf_dim);
	}

	if (y2 == 0) {
		y2 = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(y2, Buf_dim);
	}

	//------------------------------------------------

	if (error1 == 0) {
		error1 = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(error1, Buf_dim);
	}

	if (error2 == 0) {
		error2 = ippsMalloc_64f(Buf_dim);
		ippsZero_64f(error2, Buf_dim);
	}

	//------------------------------------------------

	H = new double* [M];
	for (int i = 0; i < M; i++)
	{
		H[i] = new double[N];
		memset(H[i], 0.0, (N) * sizeof(double));
	}

	F = new double* [M];
	for (int i = 0; i < M ; i++)
	{
		F[i] = new double[N];
		memset(F[i], 0.0, (N) * sizeof(double));
	}

	Hp = new double* [2*M -1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Hp[i] = new double[2*N -1];
		memset(Hp[i], 0.0, (2*N-1) * sizeof(double));
	}

	//----------------------------------------------------------------
	W1 = new double* [M];		
	for (int i = 0; i < M; i++)
	{
		W1[i] = new double[K];
		memset(W1[i], 0.0, (K) * sizeof(double));
	}

	W2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		W2[i] = new double[K];
		memset(W2[i], 0.0, (K) * sizeof(double));
	}

	W3 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		W3[i] = new double[K];
		memset(W3[i], 0.0, (K) * sizeof(double));
	}

	W4 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		W4[i] = new double[K];
		memset(W4[i], 0.0, (K) * sizeof(double));
	}
	//--------------------------------------------------
	X1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		X1[i] = new double[Buf_dim/M];
		memset(X1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	X2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		X2[i] = new double[Buf_dim / M];
		memset(X2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}
	//-------------------------------------------------------
	an_buffer1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		an_buffer1[i] = new double[N];
		memset(an_buffer1[i], 0.0, (N) * sizeof(double));
	}

	an_buffer2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		an_buffer2[i] = new double[N];
		memset(an_buffer2[i], 0.0, (N) * sizeof(double));
	}

	//-------------------------------------------------------

	h11xl_buffer = new double* [2*M-1];
	for (int i = 0; i < 2*M-1; i++)
	{
		h11xl_buffer[i] = new double[2*N-1];
		memset(h11xl_buffer[i], 0.0, (2*N-1) * sizeof(double));
	}

	h12xl_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h12xl_buffer[i] = new double[2 * N - 1];
		memset(h12xl_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	h21xl_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h21xl_buffer[i] = new double[2 * N - 1];
		memset(h21xl_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	h22xl_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h22xl_buffer[i] = new double[2 * N - 1];
		memset(h22xl_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	//---------------------------------------------------------------

	h11xr_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h11xr_buffer[i] = new double[2 * N - 1];
		memset(h11xr_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	h12xr_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h12xr_buffer[i] = new double[2 * N - 1];
		memset(h12xr_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	h21xr_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h21xr_buffer[i] = new double[2 * N - 1];
		memset(h21xr_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	h22xr_buffer = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		h22xr_buffer[i] = new double[2 * N - 1];
		memset(h22xr_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}

	//---------------------------------------------------------------

	out_buf_dly1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_buf_dly1[i] = new double[Buf_dim/M];
		memset(out_buf_dly1[i], 0.0, (Buf_dim/M) * sizeof(double));
	}
	
	out_buf_dly2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_buf_dly2[i] = new double[Buf_dim / M];
		memset(out_buf_dly2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}
	dly1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		dly1[i] = new double[Buf_dim / M];
		memset(dly1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	dly2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		dly2[i] = new double[Buf_dim / M];
		memset(dly2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}
	//----------------------------------------------------------------

	Z_w1_1 = new double* [2*M-1];
	for (int i = 0; i < 2*M-1; i++)
	{
		Z_w1_1[i] = new double[K];
		memset(Z_w1_1[i], 0.0, (K) * sizeof(double));
	}

	Z_w2_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w2_1[i] = new double[K];
		memset(Z_w2_1[i], 0.0, (K) * sizeof(double));
	}

	Z_w3_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w3_1[i] = new double[K];
		memset(Z_w3_1[i], 0.0, (K) * sizeof(double));
	}

	Z_w4_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w4_1[i] = new double[K];
		memset(Z_w4_1[i], 0.0, (K) * sizeof(double));
	}

	//----------------------------------------------------------------

	Z_w1_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w1_2[i] = new double[K];
		memset(Z_w1_2[i], 0.0, (K) * sizeof(double));
	}

	Z_w2_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w2_2[i] = new double[K];
		memset(Z_w2_2[i], 0.0, (K) * sizeof(double));
	}

	Z_w3_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w3_2[i] = new double[K];
		memset(Z_w3_2[i], 0.0, (K) * sizeof(double));
	}

	Z_w4_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Z_w4_2[i] = new double[K];
		memset(Z_w4_2[i], 0.0, (K) * sizeof(double));
	}

	//----------------------------------------------------------------

	out_w1_1 = new double* [M];
	for (int i = 0; i <M; i++)
	{
		out_w1_1[i] = new double[Buf_dim/M];
		memset(out_w1_1[i], 0.0, (Buf_dim/M) * sizeof(double));
	}

	out_w2_1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w2_1[i] = new double[Buf_dim / M];
		memset(out_w2_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_w3_1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w3_1[i] = new double[Buf_dim / M];
		memset(out_w3_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_w4_1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w4_1[i] = new double[Buf_dim / M];
		memset(out_w4_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	//----------------------------------------------------------------

	out_w1_2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w1_2[i] = new double[Buf_dim / M];
		memset(out_w1_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_w2_2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w2_2[i] = new double[Buf_dim / M];
		memset(out_w2_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_w3_2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w3_2[i] = new double[Buf_dim / M];
		memset(out_w3_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_w4_2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_w4_2[i] = new double[Buf_dim / M];
		memset(out_w4_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	//----------------------------------------------------------------

	out_sum_1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_sum_1[i] = new double[Buf_dim / M];
		memset(out_sum_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_sum_2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		out_sum_2[i] = new double[Buf_dim / M];
		memset(out_sum_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	//-----------------------------------------

	E1 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		E1[i] = new double[Buf_dim / M];
		memset(E1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	E2 = new double* [M];
	for (int i = 0; i < M; i++)
	{
		E2[i] = new double[Buf_dim / M];
		memset(E2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	//-----------------------------------------

	y1_buf = new double* [M];
	for (int i = 0; i < M; i++)
	{
		y1_buf[i] = new double[N];
		memset(y1_buf[i], 0.0, (N) * sizeof(double));
	}

	y2_buf = new double* [M];
	for (int i = 0; i < M; i++)
	{
		y2_buf[i] = new double[N];
		memset(y2_buf[i], 0.0, (N) * sizeof(double));
	}

	//----------------------------------------

	error1_buf = new double* [M];
	for (int i = 0; i < M; i++)
	{
		error1_buf[i] = new double[N];
		memset(error1_buf[i], 0.0, (N) * sizeof(double));
	}

	error2_buf = new double* [M];
	for (int i = 0; i < M; i++)
	{
		error2_buf[i] = new double[N];
		memset(error2_buf[i], 0.0, (N) * sizeof(double));
	}

	//----------------------------------------

	out_M1_1 = new double* [2*M-1];
	for (int i = 0; i < 2*M-1; i++)
	{
		out_M1_1[i] = new double[Buf_dim / M];
		memset(out_M1_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M2_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M2_1[i] = new double[Buf_dim / M];
		memset(out_M2_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M3_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M3_1[i] = new double[Buf_dim / M];
		memset(out_M3_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M4_1 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M4_1[i] = new double[Buf_dim / M];
		memset(out_M4_1[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M1_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M1_2[i] = new double[Buf_dim / M];
		memset(out_M1_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M2_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M2_2[i] = new double[Buf_dim / M];
		memset(out_M2_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M3_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M3_2[i] = new double[Buf_dim / M];
		memset(out_M3_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	out_M4_2 = new double* [2 * M - 1];
	for (int i = 0; i < 2 * M - 1; i++)
	{
		out_M4_2[i] = new double[Buf_dim / M];
		memset(out_M4_2[i], 0.0, (Buf_dim / M) * sizeof(double));
	}

	//----------------------------------------

	read_dat(save_name, p0, N);  // funzione per l'importazione del filtro prototipo

	
	cos_h(p0, H,  M,  N);       // funzione che genera la matrice di filtri di analisi
	cos_p(p0, F, M,  N);        // funzione che genera la matrice di filtri di sintesi
	petr_cos_h(p0, Hp, M, N);   // funzione che genera la matrice di filtri del banco petraglia
	//voglio solo parte reale

	//importo HRTF

	read_dat(hll_name, hll, L);
	read_dat(hlr_name, hlr, L);
	read_dat(hrl_name, hrl, L);
	read_dat(hrr_name, hrr, L);

	//normalizzo HRTF
	ippsCopy_64f(hll, hrtf, L);
	ippsCopy_64f(hlr, hrtf + 1024, L);
	ippsCopy_64f(hrl, hrtf + 2 * L, L);
	ippsCopy_64f(hrr, hrtf + 3 * L, L);
	ippsMax_64f(hrtf, 4 * L, hrtfmax);

	for (size_t i = 0; i < L; i++)
	{
		hll[i] = hll[i] / hrtfmax[0];
		hlr[i] = hlr[i] / hrtfmax[0];
		hrl[i] = hrl[i] / hrtfmax[0];
		hrr[i] = hrr[i] / hrtfmax[0];
	}


}

void __stdcall PlugIn::LEPlugin_Delete()
{


	//calcG(G, F, M, K, N);

	for (int i = 0; i < M; i++)
		delete[] H[i];
	delete[] H;

	for (int i = 0; i < M; i++)
		delete[] F[i];
	delete[] F;


	for (int i = 0; i < 2*M - 1; i++)
		delete[] Hp[i];
	delete[] Hp;

//-----------------------------------------------
	for (int i = 0; i < M; i++)
		delete[] W1[i];
	delete[] W1;

	for (int i = 0; i < M; i++)
		delete[] W2[i];
	delete[] W2;

	for (int i = 0; i < M; i++)
		delete[] W3[i];
	delete[] W3;

	for (int i = 0; i < M; i++)
		delete[] W4[i];
	delete[] W4;
//------------------------------------------
	for (int i = 0; i < M; i++)
		delete[] X1[i];
	delete[] X1;

	for (int i = 0; i < M; i++)
		delete[] X2[i];
	delete[] X2;
//------------------------------------------

	for (int i = 0; i < M; i++)
		delete[] an_buffer1[i];
	delete[] an_buffer1;

	for (int i = 0; i < M; i++)
		delete[] an_buffer2[i];
	delete[] an_buffer2;

//-----------------------------------------
	for (int i = 0; i < M; i++)
		delete[] h11xl_buffer[i];
	delete[] h11xl_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xl_buffer[i];
	delete[] h12xl_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xl_buffer[i];
	delete[] h12xl_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xl_buffer[i];
	delete[] h12xl_buffer;

	//-------------------------

	for (int i = 0; i < M; i++)
		delete[] h11xr_buffer[i];
	delete[] h11xr_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xr_buffer[i];
	delete[] h12xr_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xr_buffer[i];
	delete[] h12xr_buffer;

	for (int i = 0; i < M; i++)
		delete[] h12xr_buffer[i];
	delete[] h12xr_buffer;

	//-------------------------

	for (int i = 0; i < M; i++)
		delete[] out_buf_dly1[i];
	delete[] out_buf_dly1;

	for (int i = 0; i < M; i++)
		delete[] out_buf_dly2[i];
	delete[] out_buf_dly2;

	for (int i = 0; i < M; i++)
		delete[] dly1[i];
	delete[] dly1;

	for (int i = 0; i < M; i++)
		delete[] dly2[i];
	delete[] dly2;

	//-------------------------------

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w1_1[i];
	delete[] Z_w1_1;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w2_1[i];
	delete[] Z_w2_1;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w3_1[i];
	delete[] Z_w3_1;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w4_1[i];
	delete[] Z_w4_1;

	//-------------------------------

	for (int i = 0; i < 2*M-1; i++)
		delete[] Z_w1_2[i];
	delete[] Z_w1_2;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w2_2[i];
	delete[] Z_w2_2;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w3_2[i];
	delete[] Z_w3_2;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] Z_w4_2[i];
	delete[] Z_w4_2;

	//---------------------------

	for (int i = 0; i <M; i++)
		delete[] out_w1_1[i];
	delete[] out_w1_1;

	for (int i = 0; i < M; i++)
		delete[] out_w2_1[i];
	delete[] out_w2_1;

	for (int i = 0; i < M; i++)
		delete[] out_w3_1[i];
	delete[] out_w3_1;

	for (int i = 0; i < M; i++)
		delete[] out_w4_1[i];
	delete[] out_w4_1;


	//---------------------------

	for (int i = 0; i < M; i++)
		delete[] out_w1_2[i];
	delete[] out_w1_2;

	for (int i = 0; i < M; i++)
		delete[] out_w2_2[i];
	delete[] out_w2_2;

	for (int i = 0; i < M; i++)
		delete[] out_w3_2[i];
	delete[] out_w3_2;

	for (int i = 0; i < M; i++)
		delete[] out_w4_2[i];
	delete[] out_w4_2;

	//---------------------------

	for (int i = 0; i < M; i++)
		delete[] out_sum_1[i];
	delete[] out_sum_1;

	for (int i = 0; i < M; i++)
		delete[] out_sum_2[i];
	delete[] out_sum_2;

	//---------------------------

	for (int i = 0; i < M; i++)
		delete[] E1[i];
	delete[] E1;

	for (int i = 0; i < M; i++)
		delete[] E2[i];
	delete[] E2;

	//----------------------------

	for (int i = 0; i < M; i++)
		delete[] y1_buf[i];
	delete[] y1_buf;

	for (int i = 0; i < M; i++)
		delete[] y2_buf[i];
	delete[] y2_buf;

	//----------------------------

	for (int i = 0; i < M; i++)
		delete[] error1_buf[i];
	delete[] error1_buf;

	for (int i = 0; i < M; i++)
		delete[] error2_buf[i];
	delete[] error2_buf;

	//----------------------------

	for (int i = 0; i < M; i++)
		delete[] out_M1_1[i];
	delete[] out_M1_1;

	for (int i = 0; i < M; i++)
		delete[] out_M2_1[i];
	delete[] out_M2_1;

	for (int i = 0; i < M; i++)
		delete[] out_M3_1[i];
	delete[] out_M3_1;

	for (int i = 0; i < M; i++)
		delete[] out_M1_1[i];
	delete[] out_M1_1;

	for (int i = 0; i < M; i++)
		delete[] out_M4_1[i];
	delete[] out_M4_1;

	for (int i = 0; i < M; i++)
		delete[] out_M1_2[i];
	delete[] out_M1_2;

	for (int i = 0; i < M; i++)
		delete[] out_M2_2[i];
	delete[] out_M2_2;

	for (int i = 0; i < M; i++)
		delete[] out_M3_2[i];
	delete[] out_M3_2;

	for (int i = 0; i < M; i++)
		delete[] out_M4_2[i];
	delete[] out_M4_2;

	//----------------------------

	if (p0 != 0)
	{
		ippsFree(p0);
		p0 = 0;
	}

	if (hll != 0)
	{
		ippsFree(hll);
		hll = 0;
	}

	if (hlr != 0)
	{
		ippsFree(hlr);
		hlr = 0;
	}

	if (hrl != 0)
	{
		ippsFree(hrl);
		hrl = 0;
	}

	if (hrr != 0)
	{
		ippsFree(hrr);
		hrr = 0;
	}

	if (hrtf != 0)
	{
		ippsFree(hrtf);
		hrtf = 0;
	}

	if (hrtfmax != 0)
	{
		ippsFree(hrtfmax);
		hrtfmax = 0;
	}
	//--------------------------------

	if (h11xl != 0)
	{
		ippsFree(h11xl);
		h11xl = 0;
	}

	if (h12xl != 0)
	{
		ippsFree(h12xl);
		h12xl = 0;
	}

	if (h21xl != 0)
	{
		ippsFree(h21xl);
		h21xl = 0;
	}

	if (h22xl != 0)
	{
		ippsFree(h22xl);
		h22xl = 0;
	}

	//------------------------------

	if (h11xr != 0)
	{
		ippsFree(h11xr);
		h11xr = 0;
	}

	if (h12xr != 0)
	{
		ippsFree(h12xr);
		h12xr = 0;
	}

	if (h21xr != 0)
	{
		ippsFree(h21xr);
		h21xr = 0;
	}

	if (h22xr != 0)
	{
		ippsFree(h22xr);
		h22xr = 0;
	}
	//-------------------------------

	if (h11_ds_buffer != 0)
	{
		ippsFree(h11_ds_buffer);
		h11_ds_buffer = 0;
	}

	if (h12_ds_buffer != 0)
	{
		ippsFree(h12_ds_buffer);
		h12_ds_buffer = 0;
	}

	if (h21_ds_buffer != 0)
	{
		ippsFree(h21_ds_buffer);
		h21_ds_buffer = 0;
	}

	if (h22_ds_buffer != 0)
	{
		ippsFree(h22_ds_buffer);
		h22_ds_buffer = 0;
	}

	//-------------------------------

	if (h11_sx_buffer != 0)
	{
		ippsFree(h11_sx_buffer);
		h11_sx_buffer = 0;
	}

	if (h12_sx_buffer != 0)
	{
		ippsFree(h12_sx_buffer);
		h12_sx_buffer = 0;
	}

	if (h21_sx_buffer != 0)
	{
		ippsFree(h21_sx_buffer);
		h21_sx_buffer = 0;
	}

	if (h22_sx_buffer != 0)
	{
		ippsFree(h22_sx_buffer);
		h22_sx_buffer = 0;
	}

	//------------------------------

	if (P1_1 != 0)
	{
		ippsFree(P1_1);
		P1_1 = 0;
	}

	if (P2_1 != 0)
	{
		ippsFree(P2_1);
		P2_1 = 0;
	}

	if (P3_1 != 0)
	{
		ippsFree(P3_1);
		P3_1 = 0;
	}

	if (P4_1 != 0)
	{
		ippsFree(P4_1);
		P4_1 = 0;
	}

	//------------------------------

	if (P1_2 != 0)
	{
		ippsFree(P1_2);
		P1_2 = 0;
	}

	if (P2_2 != 0)
	{
		ippsFree(P2_2);
		P2_2 = 0;
	}

	if (P3_2 != 0)
	{
		ippsFree(P3_2);
		P3_2 = 0;
	}

	if (P4_2 != 0)
	{
		ippsFree(P4_2);
		P4_2 = 0;
	}

	//-----------------------------------
	if (mu1 != 0)
	{
		ippsFree(mu1);
		mu1 = 0;
	}

	if (mu2 != 0)
	{
		ippsFree(mu2);
		mu2 = 0;
	}

	if (mu3 != 0)
	{
		ippsFree(mu3);
		mu3 = 0;
	}

	if (mu4 != 0)
	{
		ippsFree(mu4);
		mu4 = 0;
	}

	//----------------------------

	if (e1 != 0)
	{
		ippsFree(e1);
		e1 = 0;
	}

	if (e2 != 0)
	{
		ippsFree(e2);
		e2 = 0;
	}

	//----------------------------

	if (y1 != 0)
	{
		ippsFree(y1);
		y1 = 0;
	}

	if (y2 != 0)
	{
		ippsFree(y2);
		y2 = 0;
	}
	//---------------------

	if (error1 != 0)
	{
		ippsFree(error1);
		error1 = 0;
	}

	if (error2 != 0)
	{
		ippsFree(error2);
		error2 = 0;
	}
	//---------------------

}

PlugIn::~PlugIn(void)
{

}


bool __stdcall PlugIn::LEInfoIO(int index,int type, char *StrInfo)
{
	if (type == INPUT) {
		if (index == PIN_SEGNALE_L) sprintf(StrInfo, "In Left");
		if (index == PIN_SEGNALE_R) sprintf(StrInfo, "In Right");
	} 
	if (type == OUTPUT) {
		if (index == PIN_OUT_L) sprintf(StrInfo, "Out Left");
		if (index == PIN_OUT_R) sprintf(StrInfo, "Out Right");
		if (index == PIN_ERRORE_OUT_L) sprintf(StrInfo, "Out Errore left");
		if (index == PIN_ERRORE_OUT_R) sprintf(StrInfo, "Out Errore Right");
	}
	return true;
}

int __stdcall PlugIn::LESetDefPin(int index,int type, PinType *Info)
{
	if (type==OUTPUT) 
	{
		Info->DataType=PLAYBUFFER;
		Info->Exclusive=true; 
		Info->DataLen=FrameSize;
		Info->MaxDataLen=FrameSize;
		return OUTPUT; 
	}

	if (type==INPUT) 
	{
		Info->DataType=PLAYBUFFER;
		Info->Exclusive=true; 
		Info->DataLen=FrameSize;
		Info->MaxDataLen=FrameSize;
		return INPUT; 
	}

	return -1;
}

int  __stdcall PlugIn::LEConnectionRequest(int IOType,int Index,PinType *NewType)
{
	return 0;
}

LPVOID __stdcall PlugIn::LEOnNUTechMessage(int MessageType,int MessageID,WPARAM wParam,LPARAM lParam)
{
	return 0;
}


void __stdcall PlugIn::LESetName(char *Name)
{
	strncpy(Name, NUTS_NAME, MAXNAME);
}

void __stdcall PlugIn::LESetParameter(int Index,void *Data,LPVOID bBroadCastInfo)
{
	if (Index == BANDE_ID)
	{
		M = *((int*)Data);
		CBFunction(this, NUTS_UPDATERTWATCH, BANDE_ID, 0);
	}

	if (Index == LUNGHEZZA_PROTOTIPO_ID)
	{
		N = *((int*)Data);
		CBFunction(this, NUTS_UPDATERTWATCH, LUNGHEZZA_PROTOTIPO_ID, 0);
	}

	if (Index == LUNGHEZZA_INCOGNITO_ID)
	{
		L = *((int*)Data);
		CBFunction(this, NUTS_UPDATERTWATCH, LUNGHEZZA_INCOGNITO_ID, 0);
	}

	if (Index == STEPSIZE_ID)
	{
		step_size = *((double*)Data);
		CBFunction(this, NUTS_UPDATERTWATCH, STEPSIZE_ID, 0);
	}

	if (Index == PATH_ID)
	{
		strcpy(save_name, (char*)Data);
		CBFunction(this, NUTS_UPDATERTWATCH, PATH_ID, 0);
	}
	
}

int  __stdcall PlugIn::LEGetParameter(int Index,void *Data)
{
	if (Index == BANDE_ID)
	{
		*((int*)Data) = M;
	}
	

	if (Index == LUNGHEZZA_PROTOTIPO_ID)
	{
		*((int*)Data) = N;
	}
	

	if (Index == LUNGHEZZA_INCOGNITO_ID)
	{  
		*((int*)Data) = L;
	}
	

	if (Index == STEPSIZE_ID)
	{
		*((double*)Data) = step_size;
	}
	

	if (Index == PATH_ID)
	{
		strcpy((char*)Data, save_name);
	}

	return 0;
}

void __stdcall PlugIn::LESaveSetUp()
{ 
	//CBFunction(this, NUTS_WRITEFILE, 256*sizeof(char), &save_name);
}

void __stdcall PlugIn::LELoadSetUp()
{

}

void __stdcall PlugIn::LERTWatchInit()
{
	
	WatchType NewWatch1; 
	memset(&NewWatch1, 0, sizeof(WatchType)); 								
	NewWatch1.EnableWrite = true;
	NewWatch1.LenByte = sizeof(int);
	NewWatch1.TypeVar = WTC_INT;
	NewWatch1.IDVar = BANDE_ID;
	sprintf(NewWatch1.VarName, "Bande\0");
	ExtraInfoRTEdit ExEdit1;
	memset(&ExEdit1, 0, sizeof(ExtraInfoRTEdit));
	ExEdit1.TypeExtraInfo = 1;
	ExEdit1.sizeExtraInfo = sizeof(ExtraInfoRTEdit);
	ExEdit1.EnableWheel = false;
	ExEdit1.MaxValue = 512;
	ExEdit1.MinValue = 2;
	NewWatch1.ExtraInfo = &ExEdit1;
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch1);
	

	WatchType NewWatch2;
	memset(&NewWatch2, 0, sizeof(WatchType));
	NewWatch2.EnableWrite = true;
	NewWatch2.LenByte = sizeof(int);
	NewWatch2.TypeVar = WTC_INT;
	NewWatch2.IDVar = LUNGHEZZA_PROTOTIPO_ID;
	sprintf(NewWatch2.VarName, "Lunghezza filtro protoripo\0");
	ExtraInfoRTEdit ExEdit2;
	memset(&ExEdit2, 0, sizeof(ExtraInfoRTEdit));
	ExEdit2.TypeExtraInfo = 1;
	ExEdit2.sizeExtraInfo = sizeof(ExtraInfoRTEdit);
	ExEdit2.EnableWheel = false;
	ExEdit2.MaxValue = 1024;
	ExEdit2.MinValue = 2;
	NewWatch2.ExtraInfo = &ExEdit2;
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch2);


	WatchType NewWatch3;
	memset(&NewWatch3, 0, sizeof(WatchType));
	NewWatch3.EnableWrite = true;
	NewWatch3.LenByte = sizeof(int);
	NewWatch3.TypeVar = WTC_INT;
	NewWatch3.IDVar = LUNGHEZZA_INCOGNITO_ID;
	sprintf(NewWatch3.VarName, "Lunghezza HRTF\0");
	ExtraInfoRTEdit ExEdit3;
	memset(&ExEdit3, 0, sizeof(ExtraInfoRTEdit));
	ExEdit3.TypeExtraInfo = 1;
	ExEdit3.sizeExtraInfo = sizeof(ExtraInfoRTEdit);
	ExEdit3.EnableWheel = false;
	ExEdit3.MaxValue = 1024;
	ExEdit3.MinValue = 2;
	NewWatch3.ExtraInfo = &ExEdit3;
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch3);


	WatchType NewWatch4;
	memset(&NewWatch4, 0, sizeof(WatchType));
	NewWatch4.EnableWrite = true;
	NewWatch4.LenByte = sizeof(double);
	NewWatch4.TypeVar = WTC_DOUBLE;
	NewWatch4.IDVar = STEPSIZE_ID;
	sprintf(NewWatch4.VarName, "Stepsize (max value: 1/K)\0");
	ExtraInfoRTEdit ExEdit4;
	memset(&ExEdit4, 0, sizeof(ExtraInfoRTEdit));
	ExEdit4.TypeExtraInfo = 1;
	ExEdit4.sizeExtraInfo = sizeof(ExtraInfoRTEdit);
	ExEdit4.EnableWheel = false;
	ExEdit4.MaxValue = 1.0;
	ExEdit4.MinValue = 0.0;
	NewWatch4.ExtraInfo = &ExEdit4;
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch4);


	WatchType NewWatch5;
	memset(&NewWatch5, 0, sizeof(WatchType));
	NewWatch5.EnableWrite = true;
	NewWatch5.LenByte = 256 * sizeof(char);
	NewWatch5.TypeVar = WTC_LPCHAR;
	NewWatch5.IDVar = PATH_ID;
	sprintf_s(NewWatch5.VarName, MAXCARDEBUGPLUGIN, "path+nome del filtro prototipo\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch5);

	WatchType NewWatch6;
	memset(&NewWatch6, 0, sizeof(WatchType));
	NewWatch6.EnableWrite = true;
	NewWatch6.LenByte = 256 * sizeof(char);
	NewWatch6.TypeVar = WTC_LPCHAR;
	NewWatch6.IDVar = W1_PATH;
	sprintf_s(NewWatch6.VarName, MAXCARDEBUGPLUGIN, "path W1\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch6);

	WatchType NewWatch7;
	memset(&NewWatch7, 0, sizeof(WatchType));
	NewWatch7.EnableWrite = true;
	NewWatch7.LenByte = 256 * sizeof(char);
	NewWatch7.TypeVar = WTC_LPCHAR;
	NewWatch7.IDVar = W1_NAME;
	sprintf_s(NewWatch7.VarName, MAXCARDEBUGPLUGIN, "nome W1\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch7);
//----------------------------------------------------------------------
	WatchType NewWatch8;
	memset(&NewWatch8, 0, sizeof(WatchType));
	NewWatch8.EnableWrite = true;
	NewWatch8.LenByte = 256 * sizeof(char);
	NewWatch8.TypeVar = WTC_LPCHAR;
	NewWatch8.IDVar = W2_PATH;
	sprintf_s(NewWatch8.VarName, MAXCARDEBUGPLUGIN, "path W2\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch8);

	WatchType NewWatch9;
	memset(&NewWatch9, 0, sizeof(WatchType));
	NewWatch9.EnableWrite = true;
	NewWatch9.LenByte = 256 * sizeof(char);
	NewWatch9.TypeVar = WTC_LPCHAR;
	NewWatch9.IDVar = W1_NAME;
	sprintf_s(NewWatch9.VarName, MAXCARDEBUGPLUGIN, "nome W2\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch9);
//---------------------------------------------------------------------
	WatchType NewWatch10;
	memset(&NewWatch10, 0, sizeof(WatchType));
	NewWatch10.EnableWrite = true;
	NewWatch10.LenByte = 256 * sizeof(char);
	NewWatch10.TypeVar = WTC_LPCHAR;
	NewWatch10.IDVar = W3_PATH;
	sprintf_s(NewWatch10.VarName, MAXCARDEBUGPLUGIN, "path W3\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch10);

	WatchType NewWatch11;
	memset(&NewWatch11, 0, sizeof(WatchType));
	NewWatch11.EnableWrite = true;
	NewWatch11.LenByte = 256 * sizeof(char);
	NewWatch11.TypeVar = WTC_LPCHAR;
	NewWatch11.IDVar = W3_NAME;
	sprintf_s(NewWatch11.VarName, MAXCARDEBUGPLUGIN, "nome W3\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch11);
//-----------------------------------------------------------------
	WatchType NewWatch12;
	memset(&NewWatch12, 0, sizeof(WatchType));
	NewWatch12.EnableWrite = true;
	NewWatch12.LenByte = 256 * sizeof(char);
	NewWatch12.TypeVar = WTC_LPCHAR;
	NewWatch12.IDVar = W1_PATH;
	sprintf_s(NewWatch12.VarName, MAXCARDEBUGPLUGIN, "path W4\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch12);

	WatchType NewWatch13;
	memset(&NewWatch13, 0, sizeof(WatchType));
	NewWatch13.EnableWrite = true;
	NewWatch13.LenByte = 256 * sizeof(char);
	NewWatch13.TypeVar = WTC_LPCHAR;
	NewWatch13.IDVar = W1_NAME;
	sprintf_s(NewWatch13.VarName, MAXCARDEBUGPLUGIN, "nome W4\0");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch13);


}

void __stdcall PlugIn::LESampleRateChange(int NewVal,int TrigType)
{
	if(TrigType==AUDIOPROC)
	{
		if(NewVal!=SampleRate)
		{
			SampleRate = NewVal;
		}
	}

} 

void __stdcall PlugIn::LEFrameSizeChange (int NewVal,int TrigType)
{
	if(TrigType==AUDIOPROC)
	{
		if(NewVal!=FrameSize)
		{
			FrameSize = NewVal;
		}
	}
} 

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////     LOADER      ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

LEEffect * __stdcall LoadEffect(InterfaceType _CBFunction,void *PlugRef,HWND ParentDlg)
{
	return new PlugIn(_CBFunction,PlugRef,ParentDlg);
}

int __stdcall UnLoadEffect(PlugIn *effect)
{
	delete effect;
	return TRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////     GetStartUpInfo      //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

void __stdcall LENUTSDefProps(char *NameEffect,int *Width, void *data)
{
	
	strncpy(NameEffect,NUTS_NAME,MAXNAME);
	*Width=WIDTHDEF;

	if (data!=0)
	{
		StartUpNUTSProps *Info=(StartUpNUTSProps *)data;
		Info->NumInStartUp=1;
		Info->NumOutStartUp=1;
		Info->BitMaskProc = AUDIOPROC;
		Info->BitMaskDriver = OFFLINEDRIVER | DIRECTDRIVER | ASIODRIVER;
	}
}
