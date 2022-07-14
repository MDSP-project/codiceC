#include "StdAfx.h"
#include ".\plugin.h"


PlugIn::PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg): LEEffect(_CBFunction,_PlugRef,ParentDlg)
{

	LESetNumInput(2);
	LESetNumOutput(2);
	FrameSize = CBFunction(this,NUTS_GET_FS_SR,0,(LPVOID)AUDIOPROC);
	SampleRate = CBFunction(this,NUTS_GET_FS_SR,1,(LPVOID)AUDIOPROC);

	// Inizializzazione dei vari puntatori
	p0 = 0;				// Vettore dei tappi del filtro prototipo
	P = 0;				// Vettore delle potenze nelle varie sottobande
	mu = 0;				// Vettore degli step-size per ogni sottobanda
	e = 0;				// Vettore degli errori
	
	x = 0;				// Vettore dell'ingresso
	d = 0;				// Vettore del segnale di riferimento
	y = 0;				// Vettore di uscita

	// Vettori e Matrici di appoggio per le funzioni
	tempA1 = 0;	// analisi
	tempA2 = 0;	// analisi
	U1 = 0;		// analisi
	U2 = 0;		// analisi

	tempC1 = 0; // crosstalk
	tempC2 = 0; // crosstalk
	Y_tmp = 0;	// crosstalk
	
	G_adj = 0;	// adaptation

	interp = 0;	// sintesi
	tempS1 = 0;	// sintesi
	Gw = 0;		// sintesi


	memset(save_name, 0, MAX_FILE_NAME_LENGTH * sizeof(char));
	strcpy(save_name, "C:\\Users\\Gregorio\\Documents\\MATLAB\\UNIVPM\\MDSP\\SACC_Codice Matlab\\prototipoMC.dat");

	N = 128;			// Numero di tappi del filtro prototipo
	M = 16;				// Numero di sottobande
	L = 128;			// Stima della dimensione del "System Unknown" 

	//Inizializzazione dei Filtri Adattivi e del vettore delle Potenze
	K = (N + L) / M + 1;	// Numero di tappi per ogni filtro adattivo
	delay = N / M;			// Valore del ritardo
	step_size = 0.0001;  // static_cast<double>(0.1) / K; //Valore massimo dello step size per l'adattamento
	beta = 0.99;			// Peso per la stima della potenza in ogni banda
	FrameD = FrameSize / M; // Dimensione Frame decimato
	
	i = 0;					// Variabile di appoggio per saltare l'adattamento nel primo frame
}

int __stdcall PlugIn::LEPlugin_Process(PinType **Input,PinType **Output,LPVOID ExtraInfo)
{ 
	
	double* InputData_x = ((double*)Input[0]->DataBuffer);
	double* InputData_d = ((double*)Input[1]->DataBuffer);
	double* OutputData = ((double*)Output[0]->DataBuffer);
	double* OutputDataD = ((double*)Output[1]->DataBuffer); 

	// Analisi in sottobande sia del segnale di riferimento d(n) sia del segnale in ingresso x(n)
	analisi(InputData_d, InputData_x, x, d, d_buffer, x_buffer, D, X, H, Hp, M,  N,  FrameSize,tempA1, tempA2, U1, U2);

	for (int j = 0; j < FrameD; j++)
	{
		// Calcolo delle uscite e dell'errore in ogni sottobanda
		crossfilter(X, Y, X_buffer, delay_buffer, delay, e, G, D, K, M, N, FrameD, j,tempC1, Y_tmp, tempC2);
		
		if (i>0) // Saltiamo il primo frame per l'adattamento
		{
			// Normalizzazione dei valori degli step-size in base alla stima della potenza del segnale
			calculatemu(step_size, P, X, mu, M, beta, j);
			// Adattamento dei filtri G_k(z)
			adaptation(G, G_adj, mu, e, X_buffer, K, M);	
		}
		
	} 
	// Sintesi del segnale di riferimento e dell'uscita dai filtri adattivi
	sintesi(F, output_Y, Y, M, N, FrameSize, OutputData, y, interp, tempS1, Gw);
	sintesi(F, output_D, D, M, N, FrameSize, OutputDataD, y, interp, tempS1, Gw);

	i=1;

	return COMPLETED;
}

void __stdcall PlugIn::LEPlugin_Init()
{
	// Inizializzazione dei vettori necessari al funzionamento del sistema

	if (d == 0)		// Segnale di Riferimento d(n)
	{
		d = ippsMalloc_64f(FrameSize);
		ippsZero_64f(d, FrameSize);
	}

	if (x == 0)		// Segnale di ingresso x(n)
	{
		x = ippsMalloc_64f(FrameSize);
		ippsZero_64f(x, FrameSize);
	}
	
	if (y == 0)	
	{
		y = ippsMalloc_64f(FrameSize);
		ippsZero_64f(y, FrameSize);
	}
	
	
	if (p0 == 0) {		// Vettore contenente i tappi del filtro prototipo
		p0 = ippsMalloc_64f(N);			
		ippsZero_64f(p0,N);			

	}
	
	if (P == 0) {		// Vettore delle potenze in ogni sottobanda, inizializzato a 1.0
		P = ippsMalloc_64f(2 * M - 1);
		ippsSet_64f(1.0, P, 2 * M - 1);
	}
	 
	if (mu == 0) {		// Vettore degli step size per ogni banda, inizializzato a 0
		mu = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(mu, 2 * M - 1); 
	}
	
	if (e == 0)			// Vettore degli errori in ogni sottobanda
	{
		e = ippsMalloc_64f(M);
		ippsZero_64f(e, M);			
	}

	// Allocazione dei vettori di appoggio necessari alle funzioni

	if (tempA1 == 0)	// analisi
	{
		tempA1 = ippsMalloc_64f(2 * N - 1);
		ippsZero_64f(tempA1, 2 * N - 1);
	}

	if (tempA2 == 0)	// analisi
	{
		tempA2 = ippsMalloc_64f(N);
		ippsZero_64f(tempA2, N);
	}

	if (U1 == 0)		// analisi
	{
		U1 = ippsMalloc_64f(2 * M - 1);
		ippsZero_64f(U1, 2 * M - 1);
	}

	if (U2 == 0)		// analisi
	{
		U2 = ippsMalloc_64f(M);
		ippsZero_64f(U2, M);
	}

	if (tempC1 == 0)	// crosstalk
	{
		tempC1 = ippsMalloc_64f(K);
		ippsZero_64f(tempC1, K);
	}

	if (Y_tmp == 0)		// crosstalk
	{
		Y_tmp = ippsMalloc_64f(M);
		ippsZero_64f(Y_tmp, M);
	}

	if (tempC2 == 0)	// crosstalk
	{
		tempC2 = ippsMalloc_64f(delay);
		ippsZero_64f(tempC2, delay);
	}

	G_adj = new double* [M];	// adaptation
	for (int i = 0; i < M; i++)
	{
		G_adj[i] = new double[K];
		memset(G_adj[i], 0.0, (K) * sizeof(double));
	}

	interp = new double* [M];	// sintesi
	for (int i = 0; i < M; i++)
	{
		interp[i] = new double[FrameSize];
		memset(interp[i], 0.0, (FrameSize) * sizeof(double));
	}

	if (tempS1 == 0)		// Sintesi
	{
		tempS1 = ippsMalloc_64f(N);
		ippsZero_64f(tempS1, N);
	}

	if (Gw == 0)		// Sintesi
	{
		Gw = ippsMalloc_64f(M);
		ippsZero_64f(Gw, M);
	}

	// Inizializzazione delle varie Matrici necessarie al funzionamento del sistema
	H = new double* [M];				// Banco di Analisi
	for (int i = 0; i < M; i++)
	{
		H[i] = new double[N];
		memset(H[i], 0.0, (N) * sizeof(double));	
	}

	F = new double* [M];				// Banco di Sintesi
	for (int i = 0; i < M ; i++)
	{
		F[i] = new double[N];
		memset(F[i], 0.0, (N) * sizeof(double));	
	}


	Hp = new double* [2*M -1];			// Banco di Analisi della Petraglia
	for (int i = 0; i < 2 * M - 1; i++)
	{
		Hp[i] = new double[2*N -1];
		memset(Hp[i], 0.0, (2*N-1) * sizeof(double));
	}


	G = new double* [2*M-1];			// Banco dei filtri adattivi G_k(z)
	for (int i = 0; i < 2*M-1; i++)
	{
		G[i] = new double[K];
		memset(G[i], 0.0, (K) * sizeof(double));	
	}

	delay_buffer = new double* [M];		// Buffer del ritardo
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
	

	D = new double* [M];				// Segnale di riferimento in ogni banda
	for (int i = 0; i < M ; i++)
	{
		D[i] = new double[FrameD];
		memset(D[i], 0.0, (FrameD) * sizeof(double));
	}

	Y = new double* [M];				// Uscita dai filtri adattivi in  ogni banda
	for (int i = 0; i < M; i++)
	{
		Y[i] = new double[FrameD];
		memset(Y[i], 0.0, (FrameD) * sizeof(double));
	}

	d_buffer = new double* [M];			// Buffer in ingresso al banco di analisi del segnale di riferimento
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
	
	D_buffer = new double* [M];				
	for (int i = 0; i < M; i++)
	{
		D_buffer[i] = new double[N];
		memset(D_buffer[i], 0.0, (N) * sizeof(double));
	}

	output_Y = new double* [M];				// Buffer in ingresso al banco di sintesi del sistema adattivo
	for (int i = 0; i < M; i++)
	{
		output_Y[i] = new double[N];
		memset(output_Y[i], 0.0, (N) * sizeof(double));
	}


	output_D = new double* [M];				// Buffer in ingresso al banco di sintesi del riferimento
	for (int i = 0; i < M; i++)
	{
		output_D[i] = new double[N];
		memset(output_D[i], 0.0, (N) * sizeof(double));
	}


	output_E = new double* [M];				// Buffer in ingresso al banco di sintesi per l'errore
	for (int i = 0; i < M; i++)
	{
		output_E[i] = new double[N];
		memset(output_E[i], 0.0, (N) * sizeof(double));
	}
	
	// Lettura dei tappi del filtro prototipo progettato in MatLab e salvato in "prototipoMC.dat"
	read_dat(save_name, p0, N);
	// Progettazione dei banchi di analisi, analisi/Petraglia e sintesi
	petr_cos_h(p0,Hp, M,  N);		// Banco di analisi della Petraglia
	cos_h(p0, H,  M,  N);			// Banco di analisi
	cos_p(p0, F, M,  N);			// Banco di sintesi

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

	if (d != 0)
	{
		ippsFree(d);
		d = 0;
	}

	if (x != 0)
	{
		ippsFree(x);
		x = 0;
	}

	if (y != 0)
	{
		ippsFree(y);
		y = 0;
	}
	
	// Variabili di appoggio per le funzioni

	if (tempA1 != 0)		// analisi
	{
		ippsFree(tempA1);
		tempA1 = 0;
	}

	if (tempA2 != 0)		// analisi
	{
		ippsFree(tempA2);
		tempA2 = 0;
	}

	if (U1 != 0)			// analisi
	{
		ippsFree(U1);
		U1 = 0;
	}

	if (U2 != 0)			// analisi
	{
		ippsFree(U2);
		U2 = 0;
	}

	if (tempC1 != 0)		// crosstalk
	{
		ippsFree(tempC1);
		tempC1 = 0;
	}

	if (Y_tmp != 0)			// crosstalk
	{
		ippsFree(Y_tmp);
		Y_tmp = 0;
	}

	if (tempC2 != 0)			// crosstalk
	{
		ippsFree(tempC2);
		tempC2 = 0;
	}

	for (int i = 0; i < M; i++)	// adaptation
		delete[] G_adj[i];
	delete[] G_adj;


	for (int i = 0; i < M; i++) // sintesi
		delete[] interp[i];
	delete[] interp;

	if (tempS1 != 0)			// sintesi
	{
		ippsFree(tempS1);
		tempS1 = 0;
	}

	if (Gw != 0)				// sintesi
	{
		ippsFree(Gw);
		Gw = 0;
	}
}

PlugIn::~PlugIn(void)
{

}


bool __stdcall PlugIn::LEInfoIO(int index,int type, char *StrInfo)
{
	if (type == INPUT) {
		if (index == PIN_SEGNALE_IN) sprintf(StrInfo, "In Segnale");
		if (index == PIN_RIFERIMENTO_IN) sprintf(StrInfo, "In Riferimento");
	}
	if (type == OUTPUT) {
		if (index == PIN_PETRAGLIA_OUT) sprintf(StrInfo, "Out Petraglia");
		if (index == PIN_RIFERIMENTO_OUT) sprintf(StrInfo, "Out Riferimento");
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
	sprintf(NewWatch3.VarName, "Lunghezza filtro incognito\0");
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
	sprintf_s(NewWatch5.VarName, MAXCARDEBUGPLUGIN, "path del filtro prototipo");
	CBFunction(this, NUTS_ADDRTWATCH, TRUE, &NewWatch5);

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
