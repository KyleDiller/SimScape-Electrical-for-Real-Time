/*------------------------------------------------------------
*  MATE Accelerator for SimScape Electric =
*
*  Copyright Christian Dufour 2025
------------------------------------------------------------  */
//>> mex  -largeArrayDims  MateSfun.c MatrixFun.c MateUtils.c MateCalculus.c MateAdmittancePreps.c;
// woith AVX: >>mex  -largeArrayDims  COMPFLAGS='$COMPFLAGS /arch:AVX2' MateSfun.c MatrixFun.c MateUtils.c MateCalculus.c MateAdmittancePreps.c
//ssPrint_Tf("my message ..."); better

#define S_FUNCTION_NAME MateSfun
#define S_FUNCTION_LEVEL 2

#define EPS 1e-14
#define SWAP(a,b) {real_T temp=(a);(a)=(b);(b)=temp;}
#define MAX_NUMBER_SWITCHES  32

#define TRUE 1
#define FALSE 0
#define PRINT 0

#define USE_AVX 0

#ifdef USE_AVX
#include <immintrin.h>
#include "matrix.h"
#endif



#include "simstruc.h"

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "mex.h"      /* MEX-file int_Terface mechanism */
#include "mat.h"
#endif
//#include "matrix.h"
//mxDouble *mxGetPr(const mxArray *pa);




#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "MatrixFun.h"
#include "MateCalculus.h"
#include "MateUtils.h"
#include "MateAdmittancePreps.h"




/* total number of block parameters */
#define N_PAR                12


#define SIZES   	ssGetSFcnParam(S,0)		// grandeur des matrices %nb_state nb_input nb_output nb_switches nb_nodes nb_permutations nb_Itype
#define AC   		ssGetSFcnParam(S,1)		// A
#define BC        	ssGetSFcnParam(S,2) 	// B
#define CC   		ssGetSFcnParam(S,3)		// C
#define DC   		ssGetSFcnParam(S,4)		// D
#define RSW 		ssGetSFcnParam(S,5)		// Rsw
#define SWTYPE      ssGetSFcnParam(S,6)		// Switch type
#define DISC		ssGetSFcnParam(S,7)		// discretisation
#define TS   		ssGetSFcnParam(S,8)		// pas de calcul
#define X0   		ssGetSFcnParam(S,9)		// initial state
#define U0          ssGetSFcnParam(S,10)		// initial input
#define VF          ssGetSFcnParam(S,11)		// Vforward drop of switches

//MATERT.sizes,MATERT.As  ,MATERT.Bs,MATERT.Cs,MATERT.Ds,MATERT.Rsw,MATERT.SwType,MATERT.DISC,MATERT.h,MATERT.U0,MATERT.X0
//MATERT.sizes,MATERT.As  ,MATERT.Bs,MATERT.Cs,MATERT.Ds,MATERT.Rsw,MATERT.SwType,1.0,MATERT.h,MATERT.U0,MATERT.X0


/*
* mdlInitializeSizes - initialize the sizes array
*/
static void mdlInitializeSizes(SimStruct *S)
{

    int_T nb_state,nb_input,nb_output,nb_switch,nb_node,nb_perm,nb_Itype;
    


    /* Set and Check parameter count  */
    ssSetNumSFcnParams(S, N_PAR);

    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) return;

     int_T nParam = ssGetNumSFcnParams(S);  // make parameter non-tunable
     for (int_T iParam = 0; iParam < nParam; iParam++ )  ssSetSFcnParamTunable( S, iParam, 0 );





	real_T *dims=(real_T*)mxGetPr(SIZES);
	nb_state =(int_T)dims[0];
	nb_input =(int_T)dims[1];
    nb_output=(int_T)dims[2];
	nb_switch=(int_T)dims[3];
	nb_node  =(int_T)dims[4];
	nb_perm  =(int_T)dims[5];
	nb_Itype =(int_T)dims[6];
    
    mexPrintf("states:%d  inputs:%d  outputs:%d  switches:%d nodes:%d  \n", nb_state,nb_input, nb_output,nb_switch, nb_node);
     
    double mem=nb_state*nb_state+nb_state*nb_input*2+nb_state*nb_output+nb_input*nb_output;
    for (int i=0; i<nb_switch; i++) mem=mem*2;
    mem=mem*8/1e6;  // double to MB
    mexPrintf("Estimated memory used: %e MB \n", mem);

    ssSetNumContStates(    S, 0);   /* number of continuous states           */
    ssSetNumDiscStates(    S, 0);  	/* number of discrete states*/

	/* etats: [ x[n] u[n]  logic[n]   */

	if (!ssSetNumInputPorts(S, 3)) return;  // 
	  ssSetInputPortWidth(S, 0, nb_input);               
	  ssSetInputPortWidth(S, 1, nb_node);
    if (nb_switch==0) ssSetInputPortWidth(S, 2, 1);
    else ssSetInputPortWidth(S, 2, nb_switch);
    

    ssSetInputPortDirectFeedThrough(S, 0, true);
    ssSetInputPortDirectFeedThrough(S, 1, true);
    ssSetInputPortDirectFeedThrough(S, 2, true);



	//    number of outputs 
	if (!ssSetNumOutputPorts(S, 4)) return;     
    
	ssSetOutputPortWidth(S, 0, nb_node*nb_node);
    ssSetOutputPortWidth(S, 1, nb_node);
    ssSetOutputPortWidth(S, 2, nb_output);
    ssSetOutputPortWidth(S, 3, nb_switch);  // switch status
			
/* sample times */
    ssSetNumSampleTimes(   S, 1 );
    //ssSupportsMultipleExecInstances(S, true);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    

    /* Set this S-function as runtime thread-safe for multicore execution */
    //ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);
    
    /* options */
    ssSetOptions(S,SS_OPTION_EXCEPTION_FREE_CODE);
    
    

	ssSetNumDWork( S, 3);

    ssSetDWorkWidth   (S, 0, 1);
    ssSetDWorkDataType(S, 0, SS_DOUBLE);

    ssSetDWorkWidth   (S, 1, MAX_NUMBER_SWITCHES);
    ssSetDWorkDataType(S, 1, SS_DOUBLE);

    ssSetDWorkWidth   (S, 2, nb_input);
    ssSetDWorkDataType(S, 2, SS_DOUBLE);

	ssSetNumIWork( S, 9 );

	ssSetNumPWork( S, 25);   	/* number of pointer work vector elements*/

    ssSetOperatingPointCompliance( S, USE_DEFAULT_OPERATING_POINT);
    ssSetRuntimeThreadSafetyCompliance(S, RUNTIME_THREAD_SAFETY_COMPLIANCE_TRUE);






}

/*
* mdlInitializeSampleTimes - initialize the sample times array
*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
	real_T *Tsam;
	if (mxIsEmpty(TS)){
		ssSetErrorStatus(S,"1:Erreur mdlInitializeSampleTimes: Invalid Fixed Time Step");
        return;
	}
	   Tsam=(real_T*)mxGetPr(TS);
	   ssSetSampleTime(S, 0, *Tsam);
	   ssSetOffsetTime(S, 0, 0.0);
    mexPrintf("Time step= %e us\n", 1e6*(double)*Tsam);
    
}





/*
* mdlInitializeConditions - initialize the states
*/
//#define MDL_INITIALIZE_CONDITIONS
//static void mdlInitializeConditions(SimStruct *S)
// #define MDL_START                      /* Change to #undef to remove function */
// #if defined(MDL_START)
/* Function: mdlStart ==========================================================
 * Abstract:
 *      Here we cache the state (true/false) of the XDATAEVENLYSPACED parameter.
 *      We do this primarily to illustrate how to "cache" parameter values (or
 *      information which is computed from parameter values) which do not change
 *      for the duration of the simulation (or in the generated code). In this
 *      case, rather than repeated calls to mxGetPr, we save the state once.
 *      This results in a slight increase in performance.
 */
#define MDL_START
 static void mdlStart(SimStruct *S)
 {

    real_T      *Adp, *B1dp, *B2dp, *Cdp, *Ddp;
    real_T      *Yperm, *Cinj_perm, *Dinj_perm, *z_perm, *x_perm, *B2inj_perm;
    real_T *xn,*un,*sw,*swVolt,*swVoltOld,*swOld,*IhistOld,*Yold,*SrcOld,*xtmp,*temp1,*temp2, *x_temp;

    int_T nb_state,nb_input,nb_output,nb_switch,nb_node,nb_perm,nb_Itype;


    int_T *IWork=ssGetIWork(S);

    real_T           *Ac      = (real_T*)mxGetPr(AC);
    real_T           *Bc      = (real_T*)mxGetPr(BC);
    real_T           *Cc      = (real_T*)mxGetPr(CC);
    real_T           *Dc      = (real_T*)mxGetPr(DC);
    

	real_T *dims=(real_T*)mxGetPr(SIZES);
	nb_state =(int_T)dims[0];
	nb_input =(int_T)dims[1];
    nb_output=(int_T)dims[2];
	nb_switch=(int_T)dims[3];
	nb_node  =(int_T)dims[4];
	nb_perm  =(int_T)dims[5];
	nb_Itype =(int_T)dims[6];
    
    //mexPrintf("test2:%d\n", nb_input);

    real_T *Rswitch=(real_T*)mxGetPr(RSW);
    real_T *hR= (real_T*)mxGetPr(TS);
    real_T   h= (real_T)hR[0];
    real_T *discR=(real_T*)mxGetPr(DISC);
    int_T holdtype=1;   /// ok for art5, trap, be
    int_T disc=(int_T)discR[0];

    //mexPrintf("test22:%i %e %i\n", nb_input, discR[0],disc);

    // IWORK VALUES
    //
    //0:6 nb_state nb_input nb_output nb_switch nb_node nb_perm nb_Itype
    // 7 switch permutation value
    for (int_T i=0;i<7;i++){
        //IWork[i] = (int_T)dims[i];	// store various matrix sizes
        ssSetIWorkValue(S, i, (int_T)dims[i]);
    }
    ssSetIWorkValue(S, 7, 0);  
    ssSetIWorkValue(S, 8, 0);  // tics

    Adp=      (real_T*)_aligned_malloc((nb_state*nb_state*nb_perm+1)*sizeof(real_T),32); //https://stackoverflow.com/questions/9916419/why-does-assignment-to-an-element-of-an-avx-vector-wrapper-class-object-array-pr
    B1dp=     (real_T*)malloc((nb_state*nb_input*nb_perm+1)*sizeof(real_T));
    B2dp=     (real_T*)malloc((nb_state*nb_input*nb_perm+1)*sizeof(real_T));
    Cdp=      (real_T*)malloc((nb_state*nb_output*nb_perm+1)*sizeof(real_T));
    Ddp=      (real_T*)malloc((nb_input*nb_output*nb_perm+1)*sizeof(real_T));
    Yperm=    (real_T*)malloc((nb_node*nb_node*nb_perm+1)*sizeof(real_T));
    Cinj_perm=(real_T*)malloc((nb_state*nb_node*nb_perm+1)*sizeof(real_T));
    Dinj_perm=(real_T*)malloc((nb_node*(nb_input-nb_node)*nb_perm+1)*sizeof(real_T));
    z_perm=   (real_T*)malloc((nb_Itype*nb_Itype*nb_perm+1)*sizeof(real_T));
    x_perm=   (real_T*)malloc((nb_Itype*(nb_node-nb_Itype)*nb_perm+1)*sizeof(real_T));
    B2inj_perm=(real_T*)malloc((nb_state*nb_node*nb_perm+1)*sizeof(real_T));

    //ssPrintf("Y: %p %p %p %p %p %p %p %p %p %p %p\n", Adp,B1dp,B2dp,Cdp,Ddp,Yperm,Cinj_perm, Dinj_perm, z_perm, x_perm, B2inj_perm);

    //ssPrintf(" dims: %i %i %i %i %i %i %i\n",nb_state,nb_input ,nb_output,nb_switch,nb_node,nb_perm,nb_Itype);

    int ppp=7;

    for (int_T i=0;i<nb_perm;i++){
    //for (int_T i=0;i<=3;i++){
        int_T offA, offB, offC, offD;
        offA=nb_state*nb_state*i;
        offB=nb_state*nb_input*i;
        offC=nb_output*nb_state*i;
        offD=nb_output*nb_input*i;
        int_T offY, offCi, offDi, offz, offx, offB2i;
        offY =nb_node*nb_node*i;
        offCi=nb_node*nb_state*i;
        offDi=nb_node*(nb_input-nb_node)*i;
        offz =nb_Itype*nb_Itype*i;
        offx =(nb_node-nb_Itype)*nb_Itype*i;
        offB2i=nb_state*nb_node*i;

        cpmx((Adp+offA), Ac,nb_state*nb_state );
        cpmx((B1dp+offB),Bc,nb_state*nb_input );
        cpmx((Cdp+offC), Cc,nb_output*nb_state);
        cpmx((Ddp+offD), Dc,nb_output*nb_input);
        
        // if(i==ppp){
        // mexPrintf("Ac %i:\n", i);
        // printmatMEX(Ac, nb_state, nb_state);
        // mexPrintf("------------------------------------------------------:\n");
        // }


        // need to transpose input matrice because coming out of Fortran of matlab
        // well not really, if the matrices are passed as 1D vectors
        MatePermutize (      (Adp+offA), (B1dp+offB),              (Cdp+offC), (Ddp+offD), Rswitch, i,  nb_state,  nb_input, nb_output, nb_switch);
        // if(i==ppp){
        // mexPrintf("(AdpPERMU) %i:\n", i);
        // printmatMEX((Adp+offA), nb_state, nb_state);
        // mexPrintf("oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo:\n");
        // }
        

        MateDiscretize(      (Adp+offA), (B1dp+offB), (B2dp+offB), (Cdp+offC), (Ddp+offD),  h,  disc,  holdtype, nb_state,  nb_input, nb_output);

        //  if(i==ppp){
        // mexPrintf("(AdpDISCRETE) %i:\n", i);
        // printmatMEX((Adp+offA), nb_state, nb_state);
        // mexPrintf("h %e disc %i  holdtype %i:\n", h,disc,holdtype);
        // mexPrintf("ottttttttttttttttttttttttttttttttttttttttttttttttttttt:\n");
        // }   
        
        MateAdmittancePreps( (Adp+offA), (B2dp+offB),              (Cdp+offC), (Ddp+offD), (Yperm+offY), (Cinj_perm+offCi), (Dinj_perm+offDi), (z_perm+offz), (x_perm+offx),(B2inj_perm+offB2i), nb_state,  nb_input, nb_output, nb_node,nb_Itype);
        // if(i==ppp){
        // //mexPrintf("Adx %i:\n", i);
        //printmatMEX((B2inj_perm+offB2i), nb_state, nb_node);
        // }
        // if(i==ppp){
        // mexPrintf("Y %i:wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww\n", i);
        // printmatMEX((Yperm+offY), nb_node, nb_node);
        // }

    }


    //block.Dwork(1 to 10) set here!
    

    xn          =   (real_T*)_aligned_malloc( (nb_state+1)*sizeof(real_T), 32 );
    un          =   (real_T*)malloc( (nb_input+1)*sizeof(real_T) );
    sw          =   (real_T*)malloc( (nb_input+1)*sizeof(real_T) );   // extended to nb_input 16 dec 
    swVolt      =   (real_T*)malloc((nb_switch+1)*sizeof(real_T));
    swVoltOld   =   (real_T*)malloc((nb_switch+1)*sizeof(real_T));
    swOld       =   (real_T*)malloc((nb_switch+1)*sizeof(real_T));
    IhistOld    =   (real_T*)malloc((nb_node+1)*sizeof(real_T));
    Yold        =   (real_T*)malloc((nb_node*nb_node+1)*sizeof(real_T));
    SrcOld      =   (real_T*)malloc((nb_input+1)*sizeof(real_T));
    xtmp        =   (real_T*)malloc((nb_state+1)*sizeof(real_T));
    int_T tmp=nb_input;
    if (nb_output>tmp) tmp=nb_output;
    if (nb_state> tmp) tmp=nb_state;
    temp1       =   (real_T*)malloc((nb_Itype+1)*sizeof(real_T));
    //temp2       =   (real_T*)malloc(((nb_node-nb_Itype)+1)*sizeof(real_T));
    temp2       =   (real_T*)malloc(((nb_Itype)+1)*sizeof(real_T));
     //mexPrintf("temp2 size=nb_Itype (24 sept)yesssssss!!!\n");

     x_temp       =   (real_T*)_aligned_malloc( (4)*sizeof(real_T), 32 );


// init everything to 0
    for (int_T i = 0; i < nb_state;  i++) xn[i]=0.0;
    for (int_T i = 0; i < nb_input;  i++) un[i]=0.0;
    for (int_T i = 0; i < nb_switch; i++) sw[i]=0.0;
    for (int_T i = nb_switch; i < nb_input; i++) sw[i]=1.0;  // extended part with 1.
    for (int_T i = 0; i < nb_switch; i++) swVolt[i]=0.0;
    for (int_T i = 0; i < nb_switch; i++) swVoltOld[i]=0.0;
    for (int_T i = 0; i < nb_switch; i++) swOld[i]=0.0;
    for (int_T i = 0; i < nb_node;   i++) IhistOld[i]=0.0;
    for (int_T i = 0; i < nb_node*nb_node; i++) Yold[i]=0.0;
    for (int_T i = 0; i < nb_input;  i++) SrcOld[i]=0.0;
    for (int_T i = 0; i < nb_state;  i++) xtmp[i]=0.0;
    for (int_T i = 0; i < nb_Itype;  i++) temp1[i]=0.0;
    for (int_T i = 0; i < nb_Itype;  i++) temp2[i]=0.0;  //bug 30 sept (shit! encore!)
    for (int_T i = 0; i < 4;         i++) x_temp[i]=0.0;  //for AVX storage



    ssSetPWorkValue(S, 0, Adp);
    ssSetPWorkValue(S, 1, B1dp);
    ssSetPWorkValue(S, 2, B2dp);
    ssSetPWorkValue(S, 3, Cdp);
    ssSetPWorkValue(S, 4, Ddp);
    //       mexPrintf("D: %p\n", Ddp);
    //       mexPrintf("Y: %p %p %p %p %p %p\n", Yperm,Cinj_perm, Dinj_perm, z_perm, x_perm, B2inj_perm);
    // 
    //         mexPrintf("test loop:\n");
    // printmatMEX( Yperm+(nb_node*nb_node*7), nb_node, nb_node);

    ssSetPWorkValue(S, 5, Yperm);
    ssSetPWorkValue(S, 6, Cinj_perm);
    ssSetPWorkValue(S, 7, Dinj_perm);
    ssSetPWorkValue(S, 8, z_perm);
    ssSetPWorkValue(S, 9, x_perm);
    ssSetPWorkValue(S, 10, B2inj_perm);
    ssSetPWorkValue(S, 11, xn);    
    ssSetPWorkValue(S, 12, un);
    ssSetPWorkValue(S, 13, sw);
    ssSetPWorkValue(S, 14, swVolt);
    ssSetPWorkValue(S, 15, swVoltOld);
    ssSetPWorkValue(S, 16, swOld);
    ssSetPWorkValue(S, 17, IhistOld);
    ssSetPWorkValue(S, 18, Yold);
    ssSetPWorkValue(S, 19, SrcOld);
    ssSetPWorkValue(S, 20, xtmp);
    ssSetPWorkValue(S, 21, temp1);
    ssSetPWorkValue(S, 22, temp2);
    ssSetPWorkValue(S, 23, x_temp);

    real_T* Vforward;
    real_T *vfo=(real_T*)mxGetPr(VF);
    Vforward =   (real_T*)malloc(((nb_input)+1)*sizeof(real_T));
    for (int_T i = 0; i < nb_switch; i++) Vforward[i]=vfo[0];
    for (int_T i = nb_switch; i < nb_input; i++) Vforward[i]=0.0;
    ssSetPWorkValue(S, 24, Vforward);



    // //mexPrintf("test:%d\n", nb_input);
    // 
    // real_T *SwPermNumberR = (real_T*) ssGetDWork(S,0);
    // // Initialize  
    // SwPermNumberR[0] = 0;   // Sw permutation number

    

    real_T *SwType=(real_T*)mxGetPr(SWTYPE);
    real_T *swtype = (real_T*) ssGetDWork(S,1);
    for (int_T i = 0; i < nb_switch; i++) swtype[i]=SwType[i];
    for (int_T i = nb_switch; i < MAX_NUMBER_SWITCHES; i++) swtype[i]=0.0;



    real_T *Vall = (real_T*) ssGetDWork(S,2);
    for (int_T i = 0; i < nb_input; i++) Vall[i]=0.0;


     //printmatMAT(Adp+(nb_state*nb_state*0), nb_state, nb_state,"Adp1");
   
     //printmatMAT(B1dp+(nb_state*nb_state*0),nb_state,nb_input,"B1dp1");
     //printmatMAT(B2dp+(nb_state*nb_input*0),nb_state,nb_input,"B2dp1");
     //printmatMAT(Cdp+(nb_output*nb_state*0), nb_output, nb_state,"Cdp1");
     //printmatMAT(Ddp+(nb_output*nb_input*0),nb_output,nb_input,"Ddp1");
     //printmatMAT(Yperm+(nb_node*nb_node)*0,nb_node,nb_node,"Yperm1");
     //printmatMAT(Cinj_perm+(nb_node*nb_state)*0,nb_node,nb_state,"Cinj1");
     //printmatMAT(Dinj_perm+nb_node*(nb_input-nb_node)*0,nb_node,(nb_input-nb_node),"Dinj1");
     //printmatMAT(z_perm+nb_Itype*nb_Itype*0,nb_Itype,nb_Itype,"zz1");
     //  printmatMAT(x_perm+nb_Itype*(nb_node-nb_Itype)*0,nb_Itype,(nb_node-nb_Itype),"xx1");
    // printmatMAT(B2inj_perm+nb_state*nb_node*0,nb_state,nb_node,"B2inj1");

 #if USE_AVX
    __m256d veca = _mm256_setr_pd(6.0, 6.0, 6.0, 6.0);
    __m256d vecb = _mm256_setr_pd(2.0, 2.0, 2.0, 2.0);
    __m256d vecc = _mm256_setr_pd(7.0, 7.0, 7.0, 7.0);
    __m256d res_m256 = _mm256_fmadd_pd(veca, vecb, vecc);
    double res[4];
    _mm256_storeu_pd(res, res_m256);
     mexPrintf("AVX:%e %e %e %e\n", res[0],res[1],res[2],res[3]);
 #endif


}
//#endif /*  MDL_START */

/*
* mdlOutputs - compute the outputs
*/
static void mdlOutputs(SimStruct *S, int_T tid)
{

    real_T            *Y       =  ssGetOutputPortRealSignal(S,0);  // Y
    real_T			  *Ihist   =  ssGetOutputPortRealSignal(S,1);  // Ihist
    real_T			  *yout    =  ssGetOutputPortRealSignal(S,2);  // outputs
    real_T			  *sws     =  ssGetOutputPortRealSignal(S,3);  // switch status
    real_T            *x       = ssGetRealDiscStates(S);
    InputRealPtrsType uInternal= ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType Vnode    = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType swgate   = ssGetInputPortRealSignalPtrs(S,2);

    //mexPrintf("OUTPUTpause:\n");
    //mexEvalString("pause(2);");

    

    real_T *Adp         = ( real_T*)ssGetPWork(S)[0];
    real_T *B1dp        = ( real_T*)ssGetPWork(S)[1];
    real_T *B2dp        = ( real_T*)ssGetPWork(S)[2];
    real_T *Cdp         = ( real_T*)ssGetPWork(S)[3];
    real_T *Ddp         = ( real_T*)ssGetPWork(S)[4];
    real_T *Yperm       = ( real_T*)ssGetPWork(S)[5];
    real_T *Cinj_perm   = ( real_T*)ssGetPWork(S)[6];
    real_T *Dinj_perm   = ( real_T*)ssGetPWork(S)[7];
    real_T *zz          = ( real_T*)ssGetPWork(S)[8];
    real_T *xx          = ( real_T*)ssGetPWork(S)[9];
    real_T *B2inj_perm  = ( real_T*)ssGetPWork(S)[10];
    real_T *xn          = ( real_T*)ssGetPWork(S)[11];    
    real_T *un          = ( real_T*)ssGetPWork(S)[12];
    real_T *sw          = ( real_T*)ssGetPWork(S)[13];
    real_T *swVolt      = ( real_T*)ssGetPWork(S)[14];
    real_T *swVoltOld   = ( real_T*)ssGetPWork(S)[15];
    real_T *swOld       = ( real_T*)ssGetPWork(S)[16];
    real_T *IhistOld    = ( real_T*)ssGetPWork(S)[17];
    real_T *Yold        = ( real_T*)ssGetPWork(S)[18];
    real_T *SrcOld      = ( real_T*)ssGetPWork(S)[19];
    real_T *xtmp        = ( real_T*)ssGetPWork(S)[20];
    real_T *temp1       = ( real_T*)ssGetPWork(S)[21]; // max of  nb_state nb_input nb_output
    real_T *temp2       = ( real_T*)ssGetPWork(S)[22];
    real_T *x_temp      = ( real_T*)ssGetPWork(S)[23];
    real_T *Vforward    = ( real_T*)ssGetPWork(S)[24];


    int_T *dims=ssGetIWork(S);
    const int_T nb_state =(const int_T)dims[0];
    const int_T nb_input =(const int_T)dims[1];
    const int_T nb_output=(const int_T)dims[2];
    const int_T nb_switch=(const int_T)dims[3];
    const int_T nb_node  =(const int_T)dims[4];
    const int_T nb_perm  =(const int_T)dims[5];
    const int_T nb_Itype =(const int_T)dims[6];
    int_T SwPermNumber= dims[7];
    int_T tics= dims[8];

    const int_T offsetVtype =(const int_T)(dims[1]-dims[4]+dims[6]);
    const int_T offsetItype =(const int_T)(dims[1]-dims[4]);
    const int_T nb_Vtype    =(const int_T)(dims[4]-dims[6]);  


   

    real_T *swD = (real_T *)ssGetDWork(S,1);
    const int_T swtype[MAX_NUMBER_SWITCHES]={(const int_T)swD[0],  (const int_T)swD[1],  (const int_T)swD[2],  (const int_T)swD[3],  (const int_T)swD[4],  (const int_T)swD[5],  (const int_T)swD[6],  (const int_T)swD[7],
                                           (const int_T)swD[8],  (const int_T)swD[9],  (const int_T)swD[10], (const int_T)swD[11], (const int_T)swD[12], (const int_T)swD[13], (const int_T)swD[14], (const int_T)swD[15],
                                           (const int_T)swD[16], (const int_T)swD[17], (const int_T)swD[18], (const int_T)swD[19], (const int_T)swD[20], (const int_T)swD[21], (const int_T)swD[22], (const int_T)swD[23],
                                           (const int_T)swD[24], (const int_T)swD[25], (const int_T)swD[26], (const int_T)swD[27], (const int_T)swD[28], (const int_T)swD[29], (const int_T)swD[30], (const int_T)swD[31]};

    
    real_T *Vall = (real_T*) ssGetDWork(S,2);


    
//printmatMAT(xx+nb_Itype*(nb_node-nb_Itype)*7,nb_Itype,(nb_node-nb_Itype),"xx1");

// 	  % get last switch status from the last time step.
//   sw=0;
//   for j=1:nb_switch
//       if block.Dwork(3).Data(j)>0.5
//           sw=sw+2^(nb_switch-j);
//       end
//   end
// sw=sw+1;   % 0 logic to 1 logic
    real_T *SwPermNumberR = (real_T*) ssGetDWork(S,0);
    // int_T   SwPermNumber = (int_T) SwPermNumberR[0];

   
        int_T offA, offB, offC, offD,offY, offCi, offDi, offz, offx, offB2i;
        offA=(int_T)(nb_state*nb_state*SwPermNumber);
        offB=(int_T)(nb_state*nb_input*SwPermNumber);
        offC=(int_T)(nb_output*nb_state*SwPermNumber);
        offD=(int_T)(nb_output*nb_input*SwPermNumber);
        offY=(int_T)(nb_node*nb_node*SwPermNumber);
        offCi=(int_T)(nb_node*nb_state*SwPermNumber);
        offDi=(int_T)(nb_node*(nb_input-nb_node)*SwPermNumber);
        offz =(int_T)(nb_Itype*nb_Itype*SwPermNumber);
        offx =(int_T)(nb_Vtype*nb_Itype*SwPermNumber);
        offB2i=(int_T)(nb_state*nb_node*SwPermNumber);

    tics=tics+1;  // tics
            
            // % completion of state calculation with newly available V

    if (nb_Itype==0){                   // if PortType==0
        real_T tmp;
        if (nb_state>0){                // if nb_state>0.5
            for (int_T i=0; i<nb_state; i++) {
                tmp=xtmp[i];
                for (int_T j=0; j<nb_node; j++){
                    // need to precompute B(:, some_col)!!
                    tmp+=B2inj_perm[offB2i+i*nb_node+j]* *Vnode[j];        //     block.Dwork(1).Data=x_tmp_old+B2(:,end-nb_nodal_nodes+1:end)*V(Yindex);
                    //tmp+=B2dp[offB+i*nb_input+(nb_input-nb_node)+j]* *Vnode[j];
                }
                xn[i]=tmp;
            }                                                                                      
        }
        for (int_T j=0; j<nb_input-nb_node; j++)  Vall[j]       = SrcOld[j];       // // Vall=[src_old(1:nb_input_int_Ternal); V(Yindex)];
        for (int_T j=0; j<nb_node; j++) Vall[nb_input-nb_node+j]= *Vnode[j];        // 

        //printmatMAT(B2inj_perm+offB2i,nb_state, nb_node,"B2inj=");
        //printmatMAT(B2dp+offB,nb_state, nb_input,"B2dp=");
    }
    else if  (nb_Itype==nb_node){  //elseif PortType==1
        real_T tmp;
        for (int_T i=0; i<nb_node; i++) {
            tmp=0.0;
            for (int_T j=0; j<nb_node; j++){
                tmp+=Yperm[offY+i*nb_node+j] * *Vnode[j];        //         tmp=Yold*V(Yindex)+I_histold;
            }
            temp1[i]=tmp + IhistOld[i];
        }
        if (nb_state>0){                // if nb_state>0.5
            for (int_T i=0; i<nb_state; i++) {
                tmp=xtmp[i];
                for (int_T j=0; j<nb_node; j++){
                    tmp+=B2inj_perm[offB2i+i*nb_node+j]* temp1[j];          //     block.Dwork(1).Data=x_tmp_old+B2(:,nb_input_int_Ternal+1:end)*tmp;
                    //tmp+=B2dp[offB+(nb_input-nb_node)+i*nb_node+j]* temp1[j];
                }
                xn[i]=tmp;
            }
        }
        //for (int_T j=0; j<nb_input-nb_node; j++)  Vall[j]   = SrcOld[j];       //Vall=[src_old(1:nb_input_int_Ternal); tmp ; V2];
        for (int_T i=0; i<nb_node; i++) Vall[nb_input-nb_node+i] = temp1[i];        //Vall=[src_old(1:nb_input_int_Ternal); tmp];  
        //  mexPrintf("CCCCCCCC start\n");
        //  printmatMAT(Vall,1,6,"Vall");
        // printmatMAT(temp1,1,6,"temp1");
        // printmatMAT(xtmp,1,6,"xtmp");
        //  printmatMAT(xn,1,6,"Xn");
        //  printmatMAT(IhistOld,1,6,"Ihistold");
         
    }
    else{  //mixed type
        real_T tmp;
            // V1=V(Yindex);                           Ityp Vtyp
            // V2=V1(end-NumPortTypeVolt+1:end);   Vx=[ V3   V2 ]
            // V3=V1(1:end-NumPortTypeVolt);
            // tmp=zz*(xx*V2-V3)+I_histold(1:end-NumPortTypeVolt);

            for (int_T i=0; i<nb_Itype; i++) {
                tmp=0.0;
                for (int_T j=0; j<nb_Vtype; j++){                                
                    tmp+=xx[offx+i*nb_Vtype + j]* *Vnode[nb_Itype+j];          //  xx*V2 
                }

                temp1[i] = tmp - *Vnode[i];             //temp1=(xx*V2-V3)  len(temp1)=nb_Itype
            }

            for (int_T i=0; i<nb_Itype; i++) {              //     print_Tmat((z_perm+offz), NITYPE, NITYPE);
                tmp=0.0;
                for (int_T j=0; j<nb_Itype; j++){                                
                    tmp+=zz[offz+i*nb_Itype + j]* temp1[j];          //zz*(xx*V2-V3)
                }
                temp2[i] = tmp + IhistOld[i];    //tmp=zz*(xx*V2-V3)+I_histold(1:end-NumPortTypeVolt);       
            }


           //            Yperm
//      <----NNODE-------->
//         Ityp      Vtyp
//      <-NITYPE->< ----  ><NNODE-NITYPE
//      |---------|--------| _
//      |     zz  |    xx  | NITYPE
//      |------------------| -
//      |         |        |
//      |   x2    |    y   | NNODE-NITYPE
//      |---------|--------| -
        //y=Yperm{i}(k,k); 
            
       if (nb_state>0) {                 // if nb_state>0.5
            for (int_T i=0; i<nb_state; i++) {
                tmp=xtmp[i];
                for (int_T j=0; j<nb_Itype; j++){
                    tmp+=B2inj_perm[offB2i+i*nb_node+j]* temp2[j];       // block.Dwork(1).Data=x_tmp_old+B2(:,end-nb_nodal_nodes+1:end)*[tmp;V2];
                    //tmp+=B2dp[offB+i*nb_input+(nb_input-nb_node)+j]* temp1[j];       
                }
                 for (int_T j=nb_Itype; j<nb_node; j++){
                    tmp+=B2inj_perm[offB2i+i*nb_node+j]* *Vnode[j];    
                    //tmp+=B2dp[offB+i*nb_input+(nb_input-nb_node)+j]* *Vnode[j];  // note B2inj size= nb_state x nb_node
                }
                xn[i]=tmp; 
            }                                                                                     
       }
        // displaced 23 sept...
        for (int_T j=0; j<nb_input-nb_node; j++)  Vall[j]            = SrcOld[j];       //Vall=[src_old(1:nb_input_int_Ternal); tmp ; V2];
        for (int_T j=0; j<nb_Itype; j++)          Vall[offsetItype+j]=  temp2[j];        //% concatenation of input [internal nodal]
        for (int_T j=nb_Itype; j<nb_node; j++)    Vall[offsetItype+j]= *Vnode[j];       //Vall=[src_old(1:nb_input_int_Ternal); tmp ; V2];
    }
    
    //% calculation of outputs (if no states then state loop does not enter
        for (int_T i=0; i<nb_output; i++) {
            real_T tmp=0;
            for (int_T j=0; j<nb_state; j++){
                tmp+=Cdp[offC+i*nb_state+j]* xn[j];       // //Yout=C*block.Dwork(1).Data+D*Vall;
            }
            for (int_T j=0; j<nb_input; j++){
                tmp+=Ddp[offD+i*nb_input+j]* Vall[j];
            }
            yout[i]=tmp;
        }
    

    // %block.Dwork(8).Data=Yout;
    // block.OutputPort(3).Data=Yout;  % output Yout;
    // 
    // 
    // if nb_switch>0.5
    //     block.Dwork(5).Data=block.Dwork(4).Data;
    //     block.Dwork(4).Data=Yout(1:nb_switch);
    // end
    for (int_T i=0; i<nb_switch; i++){
           swVoltOld[i] = swVolt[i];
           swVolt[i]    = yout[i];
           swOld[i]=sw[i];  //block.Dwork(6).Data=block.Dwork(3).Data; save old switch status
    }


   //  %% step 2 :compute the next switch status, next Y, next Ihist
  
    //UIn_group_ordered=Uin(MATE.UinOrder);   % source 
    //block.Dwork(9).Data=UIn_group_ordered;
    for (int_T i=0; i<(nb_input-nb_node); i++) SrcOld[i]= *uInternal[i];

     ///if nb_switch>0.5
         //SWInOrder=SWin(MATE.SwOrder);   % switch signals (C: assume they are in order for the moment)
     //end
    int_T permold;
    permold=SwPermNumber;

if (1==1){
    
    {
    //int_T pow2=1;
    unsigned int pow2 =(unsigned int)(nb_perm/2);
    int_T perm=0;
    
    
    
    for (int_T i=0; i<nb_switch; i++){
        if (swtype[i]==1){             //switch
            if ( (int_T)(*swgate[i]) == 1) sw[i]=1;
            else sw[i]=0;
        }
        else if (swtype[i]==2){       //breaker
            sw[i]=swOld[i];
            if ( (int_T)(*swgate[i]) ==1 ) sw[i]=1;
            else if (swVoltOld[i]*swVolt[i]<0) sw[i]=0;
        }
        else if (swtype[i]==3) {      //diode 
            if (swVolt[i]>=Vforward[i]) sw[i]=1; // it should be simple like this
            else sw[i]=0;
        }
        else if (swtype[i]==4) {      //thyristor
            if (swVolt[i]>=Vforward[i]  &&  (*swgate[i])>0.5 ) sw[i]=1; // it should be simple like this
            else if ( swOld[i]==1 && swVolt[i]>=Vforward[i] ) sw[i]=1;
            else sw[i]=0;

            // if (1==0){
            // if ((swVolt[i]>=(Vforward[i])) && ( ((*swgate[i])>0.5) || (sw[i]>0.5) ) )  sw[i]=1;
            // if (((2*swVolt[i]-swVoltOld[i])>=(0+Vforward[i])) && ( ((*swgate[i])>0.5) || (sw[i]>0.5) ) )  sw[i]=1;
            // //else if ((2*swVolt[i]-swVoltOld[i]) < (0+Vforward[i]) ) sw[i]=0;
            // //else if ((swVolt[i]) <= (0+Vforward[i]) ) sw[i]=0;
            // //else sw[i]=swOld[i];  
            // else sw[i]=0;
            // }
            // else if (1==1) {
            // if ((swVolt[i]>=(Vforward[i])) && ( ((*swgate[i])>0.5) || (sw[i]>0.5) ) )  sw[i]=1;
            //     //else sw[i]=0;
            // if (((2*swVolt[i]-swVoltOld[i])>=(Vforward[i])) && ( ((*swgate[i])>0.5) || (sw[i]>0.5) ) )  sw[i]=1;
            // else if  ( (swVolt[i]-Vforward[i])/0.01 < 50.0) sw[i]=0;
            // //else if ((2*swVolt[i]-swVoltOld[i]) < (0+Vforward[i]) ) sw[i]=0;
            // //else if ((swVolt[i]) < (0+Vforward[i]) ) sw[i]=0;
            // //else sw[i]=swOld[i];
            // 
            // 
            // }
            // else {
            //     //&& ((2*swVolt[i]-swVoltOld[i])>=Vforward[i]) 
            //     if (swOld[i]>0.5){
            //         if ((2*swVolt[i]-swVoltOld[i])<Vforward[i]) sw[i]=0;
            //     }
            //     else if ((swVolt[i]>=Vforward[i]) && ( *swgate[i])>0.5) sw[i]=1;
            //     else sw[i]=swOld[i];
            // }
        }
        perm+=pow2*(unsigned int)sw[i];

        //pow2=pow2*2;
        pow2=(unsigned int)(pow2/2);      // compute new permutation number
        // goto in reverse 23 sept.
    }
        SwPermNumber=(int_T)perm;
    }
    if (permold!=SwPermNumber && PRINT){
             mexPrintf("switch type: %i %i %i %i\n",swtype[0],swtype[1],swtype[2],swtype[3]);
             mexPrintf("new permutation: %i \n",SwPermNumber);
        //printmatMAT(Yperm+(nb_node*nb_node)*7,nb_node,nb_node,"Yperm1");
    }
    
}
    // if (tics>100){
    //     SwPermNumber++;
    //     if (SwPermNumber>1) SwPermNumber=0;
    //     tics=0;
    //     mexPrintf("new permutation: %i \n",SwPermNumber);
    // }
  
//    SwPermNumber=0;  tested!



 // for j=1:nb_switch
 //     switch SwType(j)
 //         case 1  %switch
 //             if SWInOrder(j)>0.5
 //                 block.Dwork(3).Data(j)=1;
 //             else
 //                 block.Dwork(3).Data(j)=0;
 //             end            
 //         case 2  % breaker
 //             if SWInOrder(j)>0.5
 //                 block.Dwork(3).Data(j)=1;
 //             else if ~isequal(sign(block.Dwork(4).Data(j)),sign(block.Dwork(5).Data(j)))
 //                       block.Dwork(3).Data(j)=0;
 //                  end
 //             end
 //      end
 //  end


// for j=1:nb_switch
//       if block.Dwork(3).Data(j)>0.5
//          %sw=sw+2^(j-1);
//          sw=sw+2^(nb_switch-j);
//       end
// end
        offA=(int_T)(nb_state*nb_state*SwPermNumber);
        offB=(int_T)(nb_state*nb_input*SwPermNumber);
        offC=(int_T)(nb_output*nb_state*SwPermNumber);
        offD=(int_T)(nb_output*nb_input*SwPermNumber);
        offY =(int_T)(nb_node*nb_node*SwPermNumber);
        offCi=(int_T)(nb_node*nb_state*SwPermNumber);
        offDi=(int_T)(nb_node*(nb_input-nb_node)*SwPermNumber);
        offz =(int_T)(nb_Itype*nb_Itype*SwPermNumber);
        offx =(int_T)((nb_node-nb_Itype)*nb_Itype*SwPermNumber);
        offB2i=(int_T)(nb_state*nb_node*SwPermNumber);

     __m256d vec_a, vec_b, vec_c, res_m256;

     //if (nb_state>0) {        no need for this condition as the loop conditions taje care fo this         // if nb_state>0.5
         
         for (int_T i=0; i<nb_state; i++) {
             real_T tmp;
             tmp=0;
             #if USE_AVX
             vec_c = _mm256_setr_pd(0.0, 0.0, 0.0, 0.0);
             {
             int_T j;
             for (j=0; j<nb_state; j=j+4) {
                 vec_a=_mm256_set_pd(Adp[offA+i*nb_state+j],Adp[offA+i*nb_state+j+1],Adp[offA+i*nb_state+j+2],Adp[offA+i*nb_state+j+3]);
                 vec_b=_mm256_set_pd(xn[j],xn[j+1],xn[j+2],xn[j+3]);
                 res_m256 = _mm256_fmadd_pd(vec_a, vec_b, vec_c);
                 vec_c= res_m256;
             }
                _mm256_storeu_pd(x_temp, res_m256);
                 tmp=x_temp[0]+x_temp[1]+x_temp[2]+x_temp[3];
         
                 //mexPrintf("CCCCCCCC j=%i\n", j);
             for (int_T jj=j-4; jj<nb_state; jj++)              tmp+= Adp[offA+i*nb_state+jj]  *xn[jj];
             }
            
             #else
             for (int_T j=0; j<nb_state; j++)                   tmp+= Adp[offA+i*nb_state+j]  *xn[j];
             #endif
             


             for (int_T j=nb_switch*0; j<nb_input; j++)         tmp+= B1dp[offB+i*nb_input+j] * Vall[j];
             //ADD UIn_group_ordered(1:nb_switch,:)=UIn_group_ordered(1:nb_switch,:).* block.Dwork(3).Data here:

             for (int_T j=0; j<nb_input-nb_node; j++)           tmp+= B2dp[offB+i*nb_input+j] * (*uInternal[j] *sw[j]);
             xtmp[i]=tmp;
         }
         for (int_T i=0; i<nb_node; i++) {
             real_T tmp;
             tmp=0;
              for (int_T j=0; j<nb_state; j++)                  tmp+=Cinj_perm[offCi+i*nb_state+j]           * xtmp[j];  
              for (int_T j=0; j<nb_input-nb_node; j++)          tmp+=Dinj_perm[offDi+i*(nb_input-nb_node)+j] * (*uInternal[j] *sw[j]);  // sw is 0 or 1 but all 1 above switch index
             IhistOld[i]=tmp;
         }

     //}
    
         // printmatMAT(xtmp,1,6,"xtmp=");
         // printmatMAT(xn,1,6,"Xn=");
         // printmatMAT(IhistOld,1,6,"Ihistold(end)=");
         // mexPrintf("CCCCCCCC end\n");
            
    // if nb_state>0.5
    // 
    //     x_tmp=A*block.Dwork(1).Data ...;
    //         +B*Vall ...;  % this is all source including nodal solution known from past history
    //         +B2(:,1:nb_input_internal)*UIn_group_ordered(1:nb_input_internal,:);  %  this is only int_Ternal source are known at the current step.
    // 
    //     I_hist=Cinj*x_tmp+Dinj*UIn_group_ordered(1:nb_input_internal,:);
    // else
    //     x_tmp=[];
    //     I_hist= Dinj*UIn_group_ordered(1:nb_input_internal,:);
    // end

    // no need to store OldY if we store the swOld! //block.Dwork(8).Data=Yvec;
    // deja calcule et mis en Dwork!                //block.Dwork(10).Data=x_tmp;
    // deja calcule et mis en Dwork                 //block.Dwork(7).Data=I_hist;
    for (int_T j=0; j<nb_input-nb_node; j++) Vall[j] = (*uInternal[j] *sw[j]);    //block.Dwork(2).Data=UIn_group_ordered;

    //block.OutputPort(1).Data=Y;  % output Y
    //block.OutputPort(2).Data=-I_hist;  % output Ihist
    for (int_T j=0; j<nb_node*nb_node; j++) Y[j]    =  Yperm[offY+j];
    for (int_T j=0; j<nb_node; j++)         Ihist[j] = 0.0-IhistOld[j];

    //sws[0]=SwPermNumber; // block.OutputPort(4).Data=block.Dwork(3).Data;  % output switch status
    for (int_T j=0; j<nb_switch; j++) sws[j] = sw[j];  ///modif 20 dec 2025
    //dims[7]=SwPermNumber;  // store new permutation

    ssSetIWorkValue(S, 8, tics);
    ssSetIWorkValue(S, 7, SwPermNumber); 
    
    

}   



static void mdlTerminate(SimStruct *S)
{

 
        // int i;
        // for (i = 0; i<ssGetNumPWork(S); i++) {
        //     if (ssGetPWorkValue(S,i) != NULL) {
        //         free(ssGetPWorkValue(S,i));
        //     }
        // }
    // {0 11} are to be _aligned_free().
        for (int i = 0; i<23; i++) {  
            if ( (i==0) || (i==11) || (i==23) ) _aligned_free(ssGetPWorkValue(S,i));
            else                                         free(ssGetPWorkValue(S,i));
        }
   

}



#ifdef	MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file int_Terface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
