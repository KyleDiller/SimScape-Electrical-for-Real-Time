/*------------------------------------------------------------
*  MATE Hbridge module
*
*  Copyright Kyle Diller 2025
------------------------------------------------------------  */


#define S_FUNCTION_NAME MateHBridgeInverter
#define S_FUNCTION_LEVEL 2

#define EPS 1e-14
#define SWAP(a,b) {real_T temp=(a);(a)=(b);(b)=temp;}
#define MAX_NUMBER_Hbridge  32

#define TRUE 1
#define FALSE 0




#include "simstruc.h"

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "mex.h"           /* MEX-file int_Terface mechanism */
#include "mat.h"
#endif
//#include "matrix.h"
//mxDouble *mxGetPr(const mxArray *pa);

//#include "MateUtils.h"


#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* total number of block parameters */
#define N_PAR                4


#define NBMODULE   	ssGetSFcnParam(S,0)		// 
#define TS   	    ssGetSFcnParam(S,1)		// 
#define RON        	ssGetSFcnParam(S,2) 	// 
#define ROFF   		ssGetSFcnParam(S,3)		// 






/*
* mdlInitializeSizes - initialize the sizes array
*/
static void mdlInitializeSizes(SimStruct *S)
{

    


    /* Set and Check parameter count  */
    ssSetNumSFcnParams(S, N_PAR);

    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

     int_T nParam = ssGetNumSFcnParams(S);  // make parameter non-tunable
     for (int_T iParam = 0; iParam < nParam; iParam++ )  ssSetSFcnParamTunable( S, iParam, 0 );





	real_T *nbm=(real_T*)mxGetPr(NBMODULE);
	const int_T nb_module =(int_T)nbm[0];

    
    //mexPrintf("states:%d  inputs:%d  outputs:%d  switches:%d nodes:%d  \n", nb_state,nb_input, nb_output,nb_switch, nb_node);
     


    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);  	/* number of discrete states*/

	/* etats: [ x[n] u[n]  logic[n]   */

	if (!ssSetNumInputPorts(S, 6)) return;  // 
        ssSetInputPortWidth(S, 0, 1);        //Vnodal
	    ssSetInputPortWidth(S, 1, nb_module);        //g1  
        ssSetInputPortWidth(S, 2, nb_module);        //g2
        ssSetInputPortWidth(S, 3, nb_module);        //g3
        ssSetInputPortWidth(S, 4, nb_module);        //g4
        ssSetInputPortWidth(S, 5, nb_module);          //Vdc of each stage


    ssSetInputPortDirectFeedThrough(S, 0, true);
    ssSetInputPortDirectFeedThrough(S, 1, true);
    ssSetInputPortDirectFeedThrough(S, 2, true);
    ssSetInputPortDirectFeedThrough(S, 3, true);
    ssSetInputPortDirectFeedThrough(S, 4, true);
    ssSetInputPortDirectFeedThrough(S, 5, true);



	//    number of outputs 
	if (!ssSetNumOutputPorts(S, 4)) return;     
    
	ssSetOutputPortWidth(S, 0, 1);  // Y
    ssSetOutputPortWidth(S, 1, 1);  //Ihist
    ssSetOutputPortWidth(S, 2, nb_module);  // Vcap
    ssSetOutputPortWidth(S, 3, nb_module);  //Idc
			
/* sample times */
    ssSetNumSampleTimes(   S, 1 );
    ssSupportsMultipleExecInstances(S, true);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    

    /* Set this S-function as runtime thread-safe for multicore execution */
    //ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    
    /* options */
    ssSetOptions(S,SS_OPTION_EXCEPTION_FREE_CODE);
    
    

	ssSetNumDWork(S, 8);

    ssSetDWorkWidth   (S, 0, nb_module);     //last g1
    ssSetDWorkDataType(S, 0, SS_DOUBLE);

    ssSetDWorkWidth   (S, 1, nb_module);     //last g2
    ssSetDWorkDataType(S, 1, SS_DOUBLE);

    ssSetDWorkWidth   (S, 2, nb_module);     //last g3
    ssSetDWorkDataType(S, 2, SS_DOUBLE);

    ssSetDWorkWidth   (S, 3, nb_module);    //last g4
    ssSetDWorkDataType(S, 3, SS_DOUBLE);

    ssSetDWorkWidth   (S, 4, 1);    //last Y
    ssSetDWorkDataType(S, 4, SS_DOUBLE);

    ssSetDWorkWidth   (S, 5, 1);    //last Ihist
    ssSetDWorkDataType(S, 5, SS_DOUBLE);

    ssSetDWorkWidth   (S, 6, 1);    //last Iout (to detect zero crossing of current output)
    ssSetDWorkDataType(S, 5, SS_DOUBLE);

    ssSetDWorkWidth   (S, 7, 1);    //High Impedance status
    ssSetDWorkDataType(S, 5, SS_DOUBLE);

	ssSetNumIWork( S, 10 );   // 

	ssSetNumPWork( S, 0);   	/* number of pointer work vector elements*/






}

/*
* mdlInitializeSampleTimes - initialize the sample times array
*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
	real_T *Tsam;
	if (mxIsEmpty(TS)){
		ssSetErrorStatus(S,"Error: Invalid Fixed Time Step");
        return;
	}
	   Tsam=(real_T*)mxGetPr(TS);
	   ssSetSampleTime(S, 0, *Tsam);
	   ssSetOffsetTime(S, 0, 0.0);
    
}





/*
* mdlInitializeConditions - initialize the states
*/
//#define MDL_INITIALIZE_CONDITIONS
//static void mdlInitializeConditions(SimStruct *S)
#define MDL_START                      /* Change to #undef to remove function */
#if defined(MDL_START)
/* Function: mdlStart ==========================================================
 * Abstract:

 */
static void mdlStart(SimStruct *S)
{

    

    real_T *nbm=(real_T*)mxGetPr(NBMODULE);
	int_T nb_module =(int_T)nbm[0];

    real_T *Tsam=(real_T*)mxGetPr(TS);

    real_T *x0 = (real_T*) ssGetDWork(S,0);
    real_T *x1 = (real_T*) ssGetDWork(S,1);
    real_T *x2 = (real_T*) ssGetDWork(S,2);
    real_T *x3 = (real_T*) ssGetDWork(S,3);
    real_T *x4 = (real_T*) ssGetDWork(S,4);
    real_T *x5 = (real_T*) ssGetDWork(S,5);
    real_T *x6 = (real_T*) ssGetDWork(S,6);
    real_T *x7 = (real_T*) ssGetDWork(S,7);
    for (int_T i = 0; i < nb_module; i++) x0[i]=0.0;
    for (int_T i = 0; i < nb_module; i++) x1[i]=0.0;
    for (int_T i = 0; i < nb_module; i++) x2[i]=0.0;
    for (int_T i = 0; i < nb_module; i++) x3[i]=0.0;
    x4[0]=0;
    x5[0]=0;
    x6[0]=0;
    x7[0]=0;


    

mexPrintf("MateHBridgeInverter: Time step= %e us  Module count: %i \n", 1e6*(double)*Tsam, nb_module);





}
#endif /*  MDL_START */

/*
* mdlOutputs - compute the outputs
*/
static void mdlOutputs(SimStruct *S, int_T tid)
{

    InputRealPtrsType Vnodal    = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType g1        = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType g2        = ssGetInputPortRealSignalPtrs(S,2);
    InputRealPtrsType g3        = ssGetInputPortRealSignalPtrs(S,3);
    InputRealPtrsType g4        = ssGetInputPortRealSignalPtrs(S,4);
    InputRealPtrsType Vdc       = ssGetInputPortRealSignalPtrs(S,5);

    real_T            *Y       =  ssGetOutputPortRealSignal(S,0);  // Y
    real_T			  *Ihist   =  ssGetOutputPortRealSignal(S,1);  // Ihist
    real_T			  *Vcapa    =  ssGetOutputPortRealSignal(S,2);  // Vcap
    real_T			  *Idc     =  ssGetOutputPortRealSignal(S,3);  // Iout


    //mexPrintf("OUTPUTpause:\n");
    //mexEvalString("pause(2);");



    real_T *g1_o     = (real_T*) ssGetDWork(S,0);
    real_T *g2_o     = (real_T*) ssGetDWork(S,1);
    real_T *g3_o     = (real_T*) ssGetDWork(S,2);
    real_T *g4_o     = (real_T*) ssGetDWork(S,3);
    real_T *Y_o      = (real_T*) ssGetDWork(S,4);
    real_T *Vhist_o  = (real_T*) ssGetDWork(S,5);
    real_T *Idc_o    = (real_T*) ssGetDWork(S,6);
    real_T *HiZ_o    = (real_T*) ssGetDWork(S,7);

    real_T *Ronn=(real_T*)mxGetPr(RON);
	const real_T Ron =Ronn[0];
    real_T *Rofff=(real_T*)mxGetPr(ROFF);
	const real_T Roff =Rofff[0];

    real_T *nbm=(real_T*)mxGetPr(NBMODULE);
    const int_T nb_module =(int_T)nbm[0];




    int_T ssx, cnt;
    int_T impulse_flag=0;
    real_T  I,Iold, Vhist, Zmmc, Ymmc, Vcapasum, Vterm, Vcapasumx;


    {
        /* Output_BEGIN */
        real_T pi=3.141592653589793;
        real_T  vsrc[3];
        real_T  R_tot_discrete[3];
        real_T  RL_discrete[3];
        int ssx,k, cnt;
        int impulse_flag=0;
        real_T  I,Iold, Vhist, Ztot, Ytot, HiZ;
        double  NegI, PosI;



        Iold=Idc_o[0]; // read old current to detect zero crossing
        HiZ=HiZ_o[0];  // hiZ status

        I=(*Vnodal[0]-Vhist_o[0])*Y_o[0];  //% I = -I_hist+Vmv*Y;   I= (Vmv-Vhist)*Ybridge
        I=(-1.0)*I; // we are using current goining out of terminal + inside this code




        // begin calculation for the current step here
        Vhist=0;  // reset Vhist
        Ztot=0;  // reset
        cnt=0;


        Ztot=2*Ron*nb_module;  // in normal active states, this is always the Total Ron value.
        for (k=0;k<nb_module;k++){
            Idc[k]=0;
            if   ((*g1[k])+(*g2[k])+(*g3[k])+(*g4[k])<0.001) {
                if (((I>0) & (Iold<0)) || ((I<0) & (Iold>0)) ) {
                    Ztot=Ztot+Roff;
                    HiZ=1;
                    
                }
                if (HiZ>0.5){
                    Ztot=Ztot+Roff;
                    HiZ=1;
                }
            }
            else HiZ=0;
            if (HiZ>0.5){
                if (*Vnodal[0]>0.0){
                    Vhist  += *Vdc[k];
                    Idc[k] = I;   // compute in case of natural rectification.
                }
                else{
                    Vhist  -= *Vdc[k];
                    Idc[k] = -I;   // compute in case of natural rectification.
                }
            }
            if (HiZ<0.5){
                NegI=(double)(I<0);
                PosI=(double)(I>0);
                Vhist  += *Vdc[k]*( *g1[k] + NegI*(1.0-(*g1[k])-(*g2[k])) );
                Vhist  -= *Vdc[k]*( *g3[k] + PosI*(1.0-(*g3[k])-(*g4[k])) );
                Idc[k] += I     *( *g1[k] + NegI*(1.0-(*g1[k])-(*g2[k])) );
                Idc[k] -= I     *( *g3[k] + PosI*(1.0-(*g3[k])-(*g4[k])) );
                HiZ=0;
            }
            // g1_o[k*4+0]=g1[k];
            // g2_o[k*4+1]=g2[k];
            // g3_o[k*4+2]=g3[k];
            // g4_o[k*4+3]=g4[k];


        }
        if (HiZ>0.5  && 1==1){
            if ( (*Vnodal[0]>Vhist)     && (*Vnodal[0]>0) ) {
                Ztot=2*Ron*nb_module;  // force low impedance
                HiZ=0;
            }
            if ( (*Vnodal[0]<Vhist) && (*Vnodal[0]<0) ){
                Ztot=2*Ron*nb_module;  // force low impedance
                HiZ=0;

            }
        }



        //mexPrintf("Hiz %e   Ztot %e Vhist %e\n", HiZ, Ztot, Vhist);
        //mexPrintf("Vdc %e %e %e %e \n", *Vdc[0],*Vdc[1],*Vdc[2],*Vdc[3]);
        Ytot=1/Ztot;
        Ihist[0]=Vhist*Ytot;  // history source
        Y[0]=Ytot;
        Vhist_o[0]=Vhist;  // save history source and Z for next step (to compute the current and update the state)
        Y_o[0]=Ytot;
        Idc_o[0]=I;
        HiZ_o[0]=HiZ;
        /* Output_END */
    }

}



static void mdlTerminate(SimStruct *S)
{


}



#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file int_Terface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
