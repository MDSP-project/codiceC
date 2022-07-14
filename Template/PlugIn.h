#pragma once
#include "LEEffect.h"
#include "ipp.h"
#include "Functions.h"


#define PLAYBUFFER_TXT "PlayBuffer"
#define	SCALAR_TXT "Scalar"
#define VECTOR_TXT "Vector"
#define VARLEN_TXT "VarLen"
			
#define INT16_TXT "Integer 16bit"
#define INT32_TXT "Integer 32bit"
#define FLT32_TXT "Float 32bit"
#define DBL64_TXT "Double 64bit"

#define DATATYPEBITMASK VECTOR|SCALAR|BUFFSTREAMING|VARLEN
#define DATALENBITMASK DATAINT32|DATAINT16|DATAINT8|DATAFLOAT32|DATADOUBLE64|CUSTOM

#define WIDTHDEF 75

#define NUTS_NAME	"HRTF"
#define BANDE_ID 6
#define LUNGHEZZA_PROTOTIPO_ID 2
#define LUNGHEZZA_INCOGNITO_ID 3
#define STEPSIZE_ID 4
#define PATH_ID 5

#define PIN_SEGNALE_IN 0
#define PIN_RIFERIMENTO_IN 1

#define PIN_PETRAGLIA_OUT 0 
#define PIN_RIFERIMENTO_OUT 1 


class PlugIn :	public LEEffect
{
public:
	PlugIn(InterfaceType _CBFunction,void * _PlugRef,HWND ParentDlg);

	~PlugIn(void);
		
	void __stdcall LESetName(char *Name);

	void __stdcall LEPlugin_Init();
	int  __stdcall LEPlugin_Process(PinType **Input,PinType **Output,LPVOID ExtraInfo);
	void __stdcall LEPlugin_Delete(void);

	bool __stdcall LEInfoIO(int index,int type, char *StrInfo);

	int __stdcall LEGetNumInput(){return LEEffect::LEGetNumInput();};  
	int __stdcall LEGetNumOutput(){return LEEffect::LEGetNumOutput();};
	//int __stdcall LESetNumInput(int Val,PinType *TypeNewIn=0){return LEGetNumInput();}; 
	//int __stdcall LESetNumOutput(int Val,PinType *TypeNewOut=0){return LEGetNumOutput();};
	int __stdcall LESetNumInput(int Val, PinType *TypeNewIn = 0) { return LEEffect::LESetNumInput(Val, TypeNewIn); };
	int __stdcall LESetNumOutput(int Val, PinType *TypeNewOut = 0) { return LEEffect::LESetNumOutput(Val, TypeNewOut); };

	void __stdcall LESetParameter(int Index,void *Data,LPVOID bBroadCastInfo);
	int  __stdcall LEGetParameter(int Index,void *Data);

	int __stdcall LESetDefPin(int index,int type, PinType *Info);

	HWND __stdcall LEGetWndSet(){return 0;};
	int __stdcall  LEWinSetStatusChange(int NewStatus){return 0;};		

	void __stdcall LESaveSetUp();
	void __stdcall LELoadSetUp();

	void __stdcall LESampleRateChange(int NewVal,int TrigType); 
	void __stdcall LEFrameSizeChange (int NewVal,int TrigType); 
	void __stdcall LEConnectionChange(int IOType,int Index,bool Connected){};
	int  __stdcall LEConnectionRequest(int IOType,int Index,PinType *NewType);
	int  __stdcall LEExtraInfoPinChange(int IOType,int Index,PinExtraInfoType ExInfo){return 0;};

	void __stdcall LERTWatchInit();
	LPVOID __stdcall LEOnNUTechMessage(int MessageType,int MessageID,WPARAM wParam,LPARAM lParam);

private:

	int FrameSize,SampleRate;
	// Vettori principali del sistema
	Ipp64f* p0, *P, *mu, *x, *d, *y, *e;
	char save_name[MAX_FILE_NAME_LENGTH];
	int N, M, L, K, delay,i;
	double** Hp, ** H, **F, **G, **delay_buffer, **D, **Y, **X;
	double** d_buffer, ** x_buffer, ** X_buffer, ** D_buffer, ** output_Y, ** output_D, ** output_E;
	double step_size, beta, FrameD;

	// Variabili di appoggio per le funzioni
	double* tempA1, * tempA2, * U1, * U2;	// analisi
	double* tempC1, * tempC2, * Y_tmp;		// crosstalk
	double** G_adj;		// adaptation
	double** interp;	// sintesi
	Ipp64f* tempS1, * Gw;	// sintesi

};