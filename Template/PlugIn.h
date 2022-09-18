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
#define PIN_ERRORE_OUT 2





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
	Ipp64f* p0, *P, *mu, *input_buffer,*y, *e; //ok
	Ipp64f* hll, *hlr, *hrl, *hrr, *hrtf, *hrtfmax; //ok
	char save_name[MAX_FILE_NAME_LENGTH], hll_name[MAX_FILE_NAME_LENGTH], hlr_name[MAX_FILE_NAME_LENGTH],hrl_name[MAX_FILE_NAME_LENGTH],hrr_name[MAX_FILE_NAME_LENGTH];
	int N, M, L, K, delay,i, Buf_dim;

	double** Hp, ** H, **F, **delay_buffer, **D, **Y, **E;

	//Definizione buffer----------------------------
	double** W1, ** W2, ** W3, ** W4;//ok
	double** X1, ** X2;//ok
	double** an_buffer1, ** an_buffer2;//ok
	//-------------------------------------------
	Ipp64f* h11xl, * h21xl, * h12xl, * h22xl; //ok
	double** h11xl_buffer, ** h21xl_buffer, ** h12xl_buffer, ** h22xl_buffer;//ok

	Ipp64f* h11xr, * h21xr, * h12xr, * h22xr; //ok
	double** h11xr_buffer, ** h21xr_buffer, ** h12xr_buffer, ** h22xr_buffer;//ok

	Ipp64f* h11_sx_buffer, * h21_sx_buffer, * h12_sx_buffer, * h22_sx_buffer;//ok
	Ipp64f* h11_ds_buffer, * h21_ds_buffer, * h12_ds_buffer, * h22_ds_buffer;//ok

	double** out_buf_dly1, ** out_buf_dly2, ** dly1, ** dly2;//ok

	double** Z_w1_1, ** Z_w2_1, ** Z_w3_1, ** Z_w4_1;//ok
	double** Z_w1_2, ** Z_w2_2, ** Z_w3_2, ** Z_w4_2;//ok

	double** out_w1_1, ** out_w2_1, ** out_w3_1, ** out_w4_1;//ok
	double** out_w1_2, ** out_w2_2, ** out_w3_2, ** out_w4_2;//ok

	Ipp64f* P1_1, * P2_1, * P3_1, * P4_1;//ok
	Ipp64f* P1_2, * P2_2, * P3_2, * P4_2;//ok

	Ipp64f* mu1, * mu2, * mu3, * mu4;//ok

	double** out_sum_1, ** out_sum_2;//ok

	Ipp64f* e1, * e2;//ok
	double** E1, ** E2;//ok



	double**X_buffer, ** D_buffer, ** output_Y, ** output_D, ** output_E;
	double step_size, beta, FrameD, alpha , delay_value;


};