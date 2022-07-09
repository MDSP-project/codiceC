#include "StdAfx.h"
#include ".\plugin.h"


PlugIn::PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg): LEEffect(_CBFunction,_PlugRef,ParentDlg)
{

	LESetNumInput(2);
	LESetNumOutput(2);
	FrameSize = CBFunction(this,NUTS_GET_FS_SR,0,(LPVOID)AUDIOPROC);
	SampleRate = CBFunction(this,NUTS_GET_FS_SR,1,(LPVOID)AUDIOPROC);	
	p0 = 0;
	P = 0;
	mu = 0;
	//x = 0;
	y = 0;
	//d = 0;
	input_buffer = 0;
	memset(save_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	strcpy(save_name, "C:\\Users\\alleg\\Desktop\\prototipoMC.dat");
	N = 128;
	M = 16;
	L = 128;
	//Inizializzazione dei Filtri Adattivi e del vettore delle Potenze
	K = (N + L) / M + 1; // Numero di tappi per ogni filtro adattivo
	delay = N / M; // Valore del ritardo
	step_size = 0.0001; //Valore massimo dello step size per l'adattamento
	beta = 0.99; // Peso per la stima della potenza in ogni banda
	FrameD = FrameSize / M; // Dimensione Frame decimato
	e = 0;
	i = 0;
}

int __stdcall PlugIn::LEPlugin_Process(PinType **Input,PinType **Output,LPVOID ExtraInfo)
{ 
	
	double* InputData_x = ((double*)Input[0]->DataBuffer);
	double* InputData_d = ((double*)Input[1]->DataBuffer);
	double* OutputData = ((double*)Output[0]->DataBuffer);
	double* OutputDataD = ((double*)Output[1]->DataBuffer); 


	analisi(InputData_d, InputData_x, d_buffer, x_buffer, D, X, H, Hp, M,  N,  FrameSize);

	for (int j = 0; j < FrameD; j++)
	{
		crossfilter(X, Y, X_buffer,delay_buffer,delay,e,G,D,K,M,N,FrameD,j);
		if (i>0)
		{
			calculatemu(step_size, P, X, mu, M, beta, j);
			adaptation(G, mu, e, X_buffer, K, M);
			
		}
		
	} 
	
	sintesi(F, output_Y, Y, M, N, FrameSize, OutputData);
	sintesi(F, output_D, D, M, N, FrameSize,OutputDataD);

	/*for (int n = 0; n <FrameSize; n++)
	{
		OutputData[n] = y[n] *32768.0 * M;
	}*/

	i=1;
	return COMPLETED;
}

void __stdcall PlugIn::LEPlugin_Init()
{
	if (p0 == 0) {
		p0 = ippsMalloc_64f(N);
		ippsZero_64f(p0,N);

	}
	
	if (P == 0) {
		P = ippsMalloc_64f(2 * M - 1);
		ippsSet_64f(1.0, P, 2 * M - 1);
	}
	 
	if (mu == 0) {
		mu = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(mu, 2 * M - 1); // Vettore degli step size per ogni banda settato a 0
	}
	
	
	if (input_buffer == 0) {
		input_buffer = ippsMalloc_64f(L);	//Buffer in ingresso al "System Unknown"
		ippsZero_64f(input_buffer,L);

	}
	
	if (e == 0)
	{
		e = ippsMalloc_64f(M);
		ippsZero_64f(e, M);
	}


	if (y == 0)
	{
		y = ippsMalloc_64f(FrameSize);
		ippsZero_64f(y, FrameSize);
	}
	
	


	read_dat(save_name, p0, N);

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


	G = new double* [2*M-1];
	for (int i = 0; i < 2*M-1; i++)
	{
		G[i] = new double[K];
		memset(G[i], 0.0, (K) * sizeof(double));
	}

	delay_buffer = new double* [M];
	for (int i = 0; i < M; i++)
	{
		delay_buffer[i] = new double[delay];
		memset(delay_buffer[i], 0.0, (delay) * sizeof(double));
	}


	X = new double* [2 * M - 1];		// Ingresso ai filtri adattivi per ogni banda
	for (int i = 0; i < 2 * M - 1; i++)
	{
		X[i] = new double[FrameD];
		memset(X[i], 0.0, (FrameD) * sizeof(double));
	}
	

	D = new double* [M];		// Segnale di riferimento in ogni banda
	for (int i = 0; i < M ; i++)
	{
		D[i] = new double[FrameD];
		memset(D[i], 0.0, (FrameD) * sizeof(double));
	}

	Y = new double* [M];		// Uscita dai filtri adattivi in  ogni banda
	for (int i = 0; i < M; i++)
	{
		Y[i] = new double[FrameD];
		memset(Y[i], 0.0, (FrameD) * sizeof(double));
	}



	d_buffer = new double* [M];		// Buffer in ingresso al banco di analisi del segnale di riferimento
	for (int i = 0; i < M; i++)
	{
		d_buffer[i] = new double[N];
		memset(d_buffer[i], 0.0, (N) * sizeof(double));
	}

	x_buffer = new double* [2 * M - 1];		// Buffer in ingresso al banco di analisi della Petraglia
	for (int i = 0; i < 2 * M - 1; i++)
	{
		x_buffer[i] = new double[2 * N - 1];
		memset(x_buffer[i], 0.0, (2 * N - 1) * sizeof(double));
	}


	X_buffer = new double* [2 * M - 1];		// Buffer in ingresso al banco di filtri adattivi
	for (int i = 0; i < 2 * M - 1; i++)
	{
		X_buffer[i] = new double[K];
		memset(X_buffer[i], 0.0, (K) * sizeof(double));
	}

	D_buffer = new double* [M];		// Buffer in ingresso al banco di filtri adattivi
	for (int i = 0; i < M; i++)
	{
		D_buffer[i] = new double[N];
		memset(D_buffer[i], 0.0, (N) * sizeof(double));
	}


	output_Y = new double* [M];		// Buffer in ingresso al banco di sintesi del sistema adattivo
	for (int i = 0; i < M; i++)
	{
		output_Y[i] = new double[N];
		memset(output_Y[i], 0.0, (N) * sizeof(double));
	}


	output_D = new double* [M];		// Buffer in ingresso al banco di sintesi del riferimento
	for (int i = 0; i < M; i++)
	{
		output_D[i] = new double[N];
		memset(output_D[i], 0.0, (N) * sizeof(double));
	}


	output_E = new double* [M];		// Buffer in ingresso al banco di sintesi per l'errore
	for (int i = 0; i < M; i++)
	{
		output_E[i] = new double[N];
		memset(output_E[i], 0.0, (N) * sizeof(double));
	}
		
	

	petr_cos_h(p0,Hp, M,  N);
	cos_h(p0, H,  M,  N);
	cos_p(p0, F, M,  N);



}

void __stdcall PlugIn::LEPlugin_Delete()
{
	for (int i = 0; i < M; i++)
		delete[] H[i];
	delete[] H;

	for (int i = 0; i < M; i++)
		delete[] F[i];
	delete[] F;

	for (int i = 0; i < M; i++)
		delete[] G[i];
	delete[] G;

	for (int i = 0; i < M; i++)
		delete[] delay_buffer[i];
	delete[] delay_buffer;

	for (int i = 0; i < 2*M - 1; i++)
		delete[] Hp[i];
	delete[] Hp;



	for (int i = 0; i < M; i++)
		delete[] D[i];
	delete[] D;

	for (int i = 0; i < M; i++)
		delete[] Y[i];
	delete[] Y;


	for (int i = 0; i < 2 * M - 1; i++)
		delete[] X[i];
	delete[] X;


	for (int i = 0; i < M; i++)
		delete[] d_buffer[i];
	delete[] d_buffer;

	for (int i = 0; i < 2*M-1; i++)
		delete[] x_buffer[i];
	delete[] x_buffer;

	for (int i = 0; i < 2 * M - 1; i++)
		delete[] X_buffer[i];
	delete[] X_buffer;


	for (int i = 0; i < M ; i++)
		delete[] D_buffer[i];
	delete[] D_buffer;


	for (int i = 0; i < M; i++)
		delete[] output_Y[i];
	delete[] output_Y;


	for (int i = 0; i < M; i++)
		delete[] output_D[i];
	delete[] output_D;

	for (int i = 0; i < M; i++)
		delete[] output_E[i];
	delete[] output_E;



	


	if (p0 != 0)
	{
		ippsFree(p0);
		p0 = 0;
	}

	if (P != 0)
	{
		ippsFree(P);
		P = 0;
	}
	if (mu != 0)
	{
		ippsFree(mu);
		mu = 0;
	}

	if (e != 0)
	{
		ippsFree(e);
		e = 0;
	}

	if (y != 0)
	{
		ippsFree(y);
		y = 0;
	}

	

	if (input_buffer != 0)
	{
		ippsFree(input_buffer);
		input_buffer = 0;
	}


	
}

PlugIn::~PlugIn(void)
{

}


bool __stdcall PlugIn::LEInfoIO(int index,int type, char *StrInfo)
{
	if(type==INPUT) sprintf(StrInfo,"In[%d]",index);
	if(type==OUTPUT) sprintf(StrInfo,"Out[%d]",index);
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

}

int  __stdcall PlugIn::LEGetParameter(int Index,void *Data)
{
	return 0;
}

void __stdcall PlugIn::LESaveSetUp()
{
	
}

void __stdcall PlugIn::LELoadSetUp()
{

}

void __stdcall PlugIn::LERTWatchInit()
{

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
