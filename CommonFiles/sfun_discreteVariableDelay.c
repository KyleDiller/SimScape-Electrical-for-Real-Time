/*
* Discrete variable time delay using interpolation.
* Times are converted to number of samples.
*
* This function implement a discrete variable time delay. The delay
* time is variable but time steps are constant.
*
*   BLOCK PARAMETERS:
*
*     - MaxDelay
*        The requested delay upper bound.
*        Trying to get a greater delay is an error.
*        Type: double; Not tunable.
*
*     - Ts
*        The discrete sample time in seconds.
*        Type: double; Not tunable.
*        Default value: 50e-6 seconds
*
*     - InitialValue
*        The value to output before attaining the delay time.
*        Type: double; Not tunable.
*        Size: One value or the width of u0
*
*     - DFT
*        Allow non-direct feedthrough for u0.
*        When unchecked and the time delay < Ts the block extrapolates
*
*   PORT PROPERTIES:
*
*   The block has two input ports and one output port.
*
*     - In
*        The signal to delay.
*        Size: Variable
*        Type: double
*
*     - D
*        The delay time.
*        Size: Variable
*        Type: double
*
*     - Out
*        The delayed signal.
*        Size: Variable
*        Type: double
*
*   BLOCK RULES:
*
*     - Width rule
*        Width of u0 and y0 are the same.
*        Width of u0 and u1 are the same or width of u1 = 1
*        Width of initial value parameter must be the same
*        as the width of "u0", or 1.
*
* Author: Olivier Tremblay , March 20011, 20-Aug-2015
* Clean-up: P.Brunelle, 01-Aug-2022
* Copyright 1997-2022 Hydro-Quebec, and The MathWorks, Inc.
*/

#define S_FUNCTION_NAME sfun_discreteVariableDelay
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#ifdef MATLAB_MEX_FILE 
     #include "mex.h"
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef TRUE
     #undef TRUE
#endif

#ifdef FALSE
     #undef FALSE
#endif

#define TRUE 1
#define FALSE 0

/* S-Function parameters */

#define NPARAMS 4

#define MAX_DELAY_IDX 0
#define MAX_DELAY(S) ssGetSFcnParam(S, MAX_DELAY_IDX)

#define SMPL_TIM_IDX 1
#define SMPL_TIM(S) ssGetSFcnParam(S, SMPL_TIM_IDX)

#define INIT_VALUE_IDX 2
#define INIT_VALUE(S) ssGetSFcnParam(S, INIT_VALUE_IDX)

#define INIT_DFT_IDX 3
#define INIT_DFT(S) ssGetSFcnParam(S, INIT_DFT_IDX)

/* Real work vector */
enum RWORK_idx_tag {
    IDX_BUFFER_SIZE,
    N_RWORKS /* Must always be last! */
};

/*Pointer work vector */
enum PWORK_idx_tag {
    PWORKS_IDX_BUFFER,
    PWORKS_IDX_WRITE_IDX,
    N_PWORKS /* Must always be last! */
};

/* define error messages */
#define ERR_INVALID_SET_INPUT_DTYPE_CALL  "Invalid call to mdlSetInputPortDataType"
#define ERR_INVALID_SET_OUTPUT_DTYPE_CALL "Invalid call to mdlSetOutputPortDataType"
#define ERR_INVALID_DTYPE "Invalid input or output port data type"
#define ERR_INVALID_PARAM "Invalid parameter. The parameter must be real, integer or boolean."
#define ERR_INVALID_DELAY_SIZE "The 'D' input size must be 1 or equal to the 'In' input size"
#define ERR_INVALID_INIT_COND_SIZE "The initial condition size must be 1 or equal to the 'In' input size"
#define ERR_INVALID_MAX_SIZE "The maximum delay must be > 0"
#define ERR_INVALID_TS "The sample time must be > 0"
#define ERR_MAX_DELAY_TOO_BIG "The maximum delay must be < 1e6 * sample time"

/* Type definitions */
typedef struct dvtdInfo_tag {
    int delaySmallerThanTsFLAG;
    int delayGreaterThanMaxDelayFLAG;
    int negativeOrInvalidDelayFLAG;
} dvtdInfo;

/* Functions */
int_T mod(int_T in, int_T size);


/* Function: mdlInitializeSizes
*
*    The sizes information is used by Simulink to determine the S-function
*    block's characteristics (number of inputs, outputs, states, etc.).
*
*    The direct feedthrough flag can be either 1=yes or 0=no. It should be
*    set to 1 if the input, "u", is used in the mdlOutput function. Setting
*    this to 0 is akin to making a promise that "u" will not be used in the
*    mdlOutput function. If you break the promise, then unpredictable results
*    will occur.
*/
static void mdlInitializeSizes(SimStruct* S) {
    int_T directFeedThrough = (int_T)*mxGetPr(INIT_DFT(S));
    ssSetNumSFcnParams(S, NPARAMS); /* Number of expected parameters */

    #if defined(MATLAB_MEX_FILE)
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        /*
        * Simulink will generate an error indicating that there is a
        * parameter mismatch if the number of expected input parameters is not equal
        * to the number of parameters entered in the dialog box
        */
        return;
    }
    #endif

    /* Parameters ARE NOT tunables */

    ssSetSFcnParamTunable(S, MAX_DELAY_IDX, FALSE);
    ssSetSFcnParamTunable(S, SMPL_TIM_IDX, FALSE);
    ssSetSFcnParamTunable(S, INIT_VALUE_IDX, FALSE);
    ssSetSFcnParamTunable(S, INIT_DFT_IDX, FALSE);

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 2)) {
        return;
    }

    ssSetInputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 0, directFeedThrough);

    ssSetInputPortWidth(S, 1, DYNAMICALLY_SIZED);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortDirectFeedThrough(S, 1, TRUE);

    if (!ssSetNumOutputPorts(S, 1)) {
        return;
    }

    ssSetOutputPortWidth(S, 0, DYNAMICALLY_SIZED);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);

    ssSetNumSampleTimes(S, 1);

    ssSetNumRWork(S, N_RWORKS);
    ssSetNumPWork(S, N_PWORKS);

    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    /* Take care when specifying exception free code - see sfuntmpl.doc */
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE));

    ssSetModelReferenceSampleTimeInheritanceRule(S, DISALLOW_SAMPLE_TIME_INHERITANCE);
    ssSetSimStateCompliance(S, USE_CUSTOM_SIM_STATE);

    ssSupportsMultipleExecInstances(S, true);
}


static void mdlInitializeSampleTimes(SimStruct* S) {
    real_T* sampleTimePt = (real_T*)mxGetPr(SMPL_TIM(S));
    ssSetSampleTime(S, 0, *sampleTimePt);
    ssSetOffsetTime(S, 0, 0.0);
}


#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_WIDTH
static void mdlSetInputPortWidth(SimStruct* S, int_T port, int_T inputPortWidth) {
    ssSetInputPortWidth(S, port, inputPortWidth);
    if (port == 0) {
        ssSetOutputPortWidth(S, port, inputPortWidth);
    }
}

#define MDL_SET_OUTPUT_PORT_WIDTH
static void mdlSetOutputPortWidth(SimStruct* S, int_T port, int_T outputPortWidth) {
    ssSetInputPortWidth(S, 0, outputPortWidth);
    ssSetInputPortWidth(S, 1, outputPortWidth);
    ssSetOutputPortWidth(S, port, outputPortWidth);
}

#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
/* Function: mdlSetDefaultPortDimensionInfo
*
* In case no ports were specified, the default is an input port of width 2
* and an output port of width 1.
*/
static void mdlSetDefaultPortDimensionInfo(SimStruct* S) {
    if (ssGetInputPortWidth(S, 0) == DYNAMICALLY_SIZED) {
        ssSetInputPortWidth(S, 0, 1);
        ssSetOutputPortWidth(S, 0, 1);
    }
    if (ssGetInputPortWidth(S, 1) == DYNAMICALLY_SIZED) {
        ssSetInputPortWidth(S, 1, 1);
    }
}
#endif

#define MDL_START
#if defined(MDL_START)

static const char* allocErrorMess = "Memory allocation error in mdlStart";
static void*       sf_alloc(SimStruct* S, size_t num_elem, size_t elem_size) {
    if (ssGetErrorStatus(S) == (char_T*)0 && num_elem > 0 && elem_size > 0) {
        void* mem = calloc(num_elem, elem_size);
        if (mem == (void*)0) {
            ssSetErrorStatus(S, allocErrorMess);
        }
        return mem;
    } else {
        return (void*)0;
    }
}


static void mdlStart(SimStruct* S) {
    int_T       i, bufSize;
    real_T**    buffer;
    int_T*      idxWrite;
    int_T       uWidth        = ssGetInputPortWidth(S, 0);
    int_T       dWidth        = ssGetInputPortWidth(S, 1);
    real_T      Ts            = *mxGetPr(SMPL_TIM(S));
    real_T      maxDelay      = *mxGetPr(MAX_DELAY(S));
    int_T       initValueNbr  = (int_T)mxGetN(INIT_VALUE(S));
    int_T       initValueSize = (int_T)mxGetM(INIT_VALUE(S));
    int_T       maxBufSize    = 1000000;
    dvtdInfo*   pinfo;

    /* Check if the delay size matches with the input size */
    if ((dWidth != uWidth) && (dWidth != 1)) {
        ssSetErrorStatus(S, ERR_INVALID_DELAY_SIZE);
        return;
    }

    /* Check if the initial values are coherent with input size */
    if ((initValueNbr != uWidth) && (initValueNbr != 1)) {
        ssSetErrorStatus(S, ERR_INVALID_INIT_COND_SIZE);
        return;
    }

    /* Check if the initial values are coherent with input size */
    if (maxDelay <= 0.0) {
        ssSetErrorStatus(S, ERR_INVALID_MAX_SIZE);
        return;
    }

    if (Ts <= 0.0) {
        ssSetErrorStatus(S, ERR_INVALID_TS);
        return;
    }

    pinfo = (dvtdInfo*)malloc(sizeof(dvtdInfo));
    if (pinfo != NULL) {
        pinfo->delaySmallerThanTsFLAG       = 0;
        pinfo->delayGreaterThanMaxDelayFLAG = 0;
        pinfo->negativeOrInvalidDelayFLAG   = 0;
    }
    ssSetUserData(S, (void*)pinfo);

    /* Determine buffer size */
    if (maxDelay > maxBufSize*Ts) {
        ssSetErrorStatus(S, ERR_MAX_DELAY_TOO_BIG);
        return;
    }
    bufSize = (int_T)ceil(maxDelay / Ts) + 1;

    if (initValueSize > bufSize) {
        bufSize = initValueSize;

        #if defined(MATLAB_MEX_FILE)
        mexWarnMsgIdAndTxt("SpecializedPowerSystems:DiscreteVaraiableTimeDelay:InvalidParameter",
                           "Block %s: Delay buffer has been increased to %d samples, in order to "
                           "respect the initial condition length.",
                           ssGetPath(S), initValueSize);
        #endif

    }
    ssSetRWorkValue(S, IDX_BUFFER_SIZE, (real_T)bufSize);

    /* Allocate the buffer memory */
    buffer = (real_T**)sf_alloc(S, uWidth, sizeof(real_T*));
    for (i = 0; i < uWidth; i++) {
        buffer[i] = (real_T*)sf_alloc(S, bufSize, sizeof(real_T));
    }
    ssGetPWork(S)[PWORKS_IDX_BUFFER] = (void*)buffer;

    /* Allocate the idxWrite memory */
    idxWrite      = (int_T*)sf_alloc(S, uWidth, sizeof(int_T));
    ssGetPWork(S)[PWORKS_IDX_WRITE_IDX] = (void*)idxWrite;
}
#endif /* MDL_START */

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions */

static void mdlInitializeConditions(SimStruct* S) {
    int_T    i, j;
    real_T** buffer   = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*   idxWrite = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    int_T    bufSize  = (int_T)ssGetRWorkValue(S, IDX_BUFFER_SIZE);

    int_T uWidth = ssGetInputPortWidth(S, 0);

    int_T   directFeedThrough = (int_T)*mxGetPr(INIT_DFT(S));
    real_T* initValues        = (real_T*)mxGetPr(INIT_VALUE(S));
    int_T   initValueNbr      = (int_T)mxGetN(INIT_VALUE(S));
    int_T   initValueSize     = (int_T)mxGetM(INIT_VALUE(S));

    for (i = 0; i < uWidth; i++) {
        memset(buffer[i], 0, bufSize * sizeof(real_T));

        if ((initValueNbr == 1) && (initValueSize == 1)) {
            for (j = 0; j < bufSize; j++) {
                buffer[i][j] = initValues[0];
            }
        } else if (initValueNbr == 1) {
            memcpy(&buffer[i][bufSize - initValueSize], &initValues[0],
                   initValueSize * sizeof(real_T));
        } else {
            memcpy(&buffer[i][bufSize - initValueSize], &initValues[i * initValueSize],
                   initValueSize * sizeof(real_T));
        }
        idxWrite[i] = 0;
    }
}
/* End of MDL_INITIALIZE_CONDITIONS */


#undef MDL_GET_TIME_OF_NEXT_VAR_HIT
#if defined(MDL_GET_TIME_OF_NEXT_VAR_HIT) && (defined(MATLAB_MEX_FILE) || defined(NRT))

/* Function: mdlGetTimeOfNextVarHit
* Abstract:
*    This function is called to get the time of the next variable sample
*    time hit. This function is called once for every major integration time
*    step. It must return time of next hit by using ssSetTNext. The time of
*    the next hit must be greater than ssGetT(S).
*
*    Note, the time of next hit can be a function of the input signal(s).
*/
static void mdlGetTimeOfNextVarHit(SimStruct* S) {
    time_T timeOfNextHit = ssGetT(S) /* + offset */;
    ssSetTNext(S, timeOfNextHit);
}
#endif /* MDL_GET_TIME_OF_NEXT_VAR_HIT */


/* Function: mdlOutputs
*  y = uBuffer[indDelay]
*/
static void mdlOutputs(SimStruct* S, int_T tid) {
    int_T             bufSize, idxRead, idxWrite, i;
    real_T            u1, u2, d, corr;
    int_T             uWidth   = ssGetInputPortWidth(S, 0);
    int_T             dWidth   = ssGetInputPortWidth(S, 1);
    real_T**          buffer   = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*            idxWrt   = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    dvtdInfo*         pinfo    = (dvtdInfo*)ssGetUserData(S);
    real_T*           output   = (real_T*)ssGetOutputPortRealSignal(S, 0);
    real_T            Ts       = *mxGetPr(SMPL_TIM(S));
    real_T            maxDelay = *mxGetPr(MAX_DELAY(S));
    int_T             directFeedThrough = (int_T)*mxGetPr(INIT_DFT(S));
    InputRealPtrsType dPtrs    = ssGetInputPortRealSignalPtrs(S, 1);
    InputRealPtrsType uPtrs    = ssGetInputPortRealSignalPtrs(S, (directFeedThrough == 1 ? 0 : 1));

    bufSize = (int_T)ssGetRWorkValue(S, IDX_BUFFER_SIZE);
    for (i = 0; i < uWidth; i++) {
        idxWrite = idxWrt[i];

        /* We can get and save the current input */
        if (directFeedThrough) {
            /* Reset the indexes of the circular buffer */
            if (idxWrite >= bufSize) {
                idxWrite = 0;
            }

            /* Write the input value */
            buffer[i][idxWrite] = *uPtrs[i];

            /* Store the next write indexes */
            idxWrt[i] = idxWrite + 1;
        }

        /* Get the delay corresponding to the Input */
        if (dWidth == uWidth) {
            d = *dPtrs[i];
        } else if (dWidth == 1) {
            d = *dPtrs[0];
        } else {
            d = 0.0;
        }

        /* Correct the possible imprecision when near a multiple of the time step */
        d = d + mxGetEps();

        if ((d > maxDelay) || mxIsInf(d)) {
            /* delayGreaterThanMaxDelayFLAG needs to be set to 1*/
            pinfo->delayGreaterThanMaxDelayFLAG = 1;
            d = maxDelay;
        } else if ((d < 0) || (d != d)) {
            /* negativeOrInvalidDelayFLAG needs to be set to 1*/
            pinfo->negativeOrInvalidDelayFLAG = 1;
            d = 0;
        }

        /*Output the current value */

        /* Get the Read index */
        idxRead = (int_T)ceil((real_T)idxWrite - d / Ts);

        /* Calculate the interpolation's correction factor */
        corr = fmod(d / Ts, 1.0);

        if (corr >= 1.0) {
            corr = 0.0;
        }

        /* extrapolation  for non direct feedthrough*/
        if ((directFeedThrough == 0) && (d < Ts)) {
            /* Set the delaySmallerThanTsFLAG to 1 so we can throw a warning in mdlTerminate */
            pinfo->delaySmallerThanTsFLAG = 1;
            u1        = buffer[i][(int_T)mod(idxRead - 1, bufSize)];
            u2        = buffer[i][(int_T)mod(idxRead - 2, bufSize)];
            output[i] = u1 - (u2 - u1) * (1.0 - corr);
        }

        /* Interpolation (normal condition) */
        else {
            u1        = buffer[i][(int_T)mod(idxRead, bufSize)];
            u2        = buffer[i][(int_T)mod(idxRead - 1, bufSize)];
            output[i] = u1 + (u2 - u1) * corr;
        }
    }
}

#define MDL_UPDATE
#if defined(MDL_UPDATE)
/* Function: mdlUpdate
*
*    Called once for every major integration time step.
*    Discrete states are typically updated here, but this function is useful
*    for performing any tasks that should only take place once per
*    integration step.
*/
static void mdlUpdate(SimStruct* S, int_T tid) {
    int_T             i, bufSize, idxWrite;
    real_T**          buffer            = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*            idxWrt            = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    int_T             uWidth            = ssGetInputPortWidth(S, 0);
    InputRealPtrsType uPtrs             = ssGetInputPortRealSignalPtrs(S, 0);
    int_T             directFeedThrough = (int_T)*mxGetPr(INIT_DFT(S));

    /* Get and save the current input */
    if (directFeedThrough == 0) {
        bufSize = (int_T)ssGetRWorkValue(S, IDX_BUFFER_SIZE);
        for (i = 0; i < uWidth; i++) {
            idxWrite = idxWrt[i];

            /* Reset the indexes of the circular buffer */
            if (idxWrite >= bufSize) {
                idxWrite = 0;
            }

            /* Write the input value */
            buffer[i][idxWrite] = *uPtrs[i];

            /* Store the next write indexes */
            idxWrt[i] = idxWrite + 1;
        }
    }
}
#endif /* MDL_UPDATE */

int_T mod(int_T in, int_T size) {
    if (in >= size) {
        return in - size;
    } else if (in < 0) {
        return in + size;
    } else {
        return in;
    }
}

/* Function: mdlTerminate
*  Called when the simulation is terminated.
*/
static void mdlTerminate(SimStruct* S) {
    int       i;
    real_T    maxDelay = *mxGetPr(MAX_DELAY(S));
    real_T**  buffer   = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*    idxWrite = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    dvtdInfo* pinfo    = (dvtdInfo*)ssGetUserData(S);
    int_T       uWidth        = ssGetInputPortWidth(S, 0);

    /* If delaySmallerThanTsFLAG is 1, throw a warning */
    if (pinfo->delaySmallerThanTsFLAG == 1) {

        #if defined(MATLAB_MEX_FILE)
        mexWarnMsgIdAndTxt("SpecializedPowerSystems:DiscreteVaraiableTimeDelay:InvalidParameter",
                           "Block %s: The specified delay time is smaller than the specified "
                           "Sample time (Ts). This may cause inaccurate simulation results.",
                           ssGetPath(S));
        #endif

    }
    /* If delayGreaterThanMaxDelayFLAG is 1, throw a warning */
    if (pinfo->delayGreaterThanMaxDelayFLAG == 1) {

        #if defined(MATLAB_MEX_FILE)
        mexWarnMsgIdAndTxt("SpecializedPowerSystems:DiscreteVaraiableTimeDelay:InvalidParameter",
                           "Block %s: The delay used is greater than the specified maximum delay "
                           "parameter (%g). Limiting delay value to maximum delay.",
                           ssGetPath(S), maxDelay);
        #endif

    }
    /* If negativeOrInvalidDelayFLAG is 1, throw a warning */
    if (pinfo->negativeOrInvalidDelayFLAG == 1) {

        #if defined(MATLAB_MEX_FILE)
        mexWarnMsgIdAndTxt(
            "SpecializedPowerSystems:DiscreteVaraiableTimeDelay:InvalidParameter",
            "Block %s: Negative/NaN time delay is not allowed. Delay has been set to zero.",
            ssGetPath(S));
        #endif

    }

    /* Release pointer work vectors */
    for (i = 0; i < uWidth; i++) {
        free(buffer[i]);
    }
    free(buffer);
    free(idxWrite);

    /* Release S-Function user data */
    free(pinfo);

} /* end mdlTerminate */


/* Functions required to support SimStates */
static const char* fieldNames[] = {
    "PWORKS_IDX_BUFFER",
    "PWORKS_IDX_WRITE_IDX",
    "IDX_BUFFER_SIZE",
};
static const int nFields = sizeof(fieldNames) / sizeof(const char*);

/* Define to indicate that this S-function has the mdlG[S]etSimState method */
#if defined(MATLAB_MEX_FILE)
#define MDL_SIM_STATE

static mxArray* mdlGetSimState(SimStruct* S) {
    int_T    i, j;
    real_T** buffer   = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*   idxWrt   = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    int_T    uWidth  = (int_T)ssGetInputPortWidth(S, 0);
    int_T    bufSize = (int_T)ssGetRWorkValue(S, IDX_BUFFER_SIZE);

    /* Create a MATLAB structure to hold the run-time data */
    mxArray* simSnap = mxCreateStructMatrix(1, 1, nFields, fieldNames);
    {
        mxArray* mPWork   = mxCreateDoubleMatrix(1, uWidth * bufSize, mxREAL);
        real_T*  mPWorkMx = mxCalloc(uWidth * bufSize, sizeof(real_T));
        for (i = 0; i < uWidth; i++) {
            for (j = 0; j < bufSize; j++) {
                mPWorkMx[i * bufSize + j] = buffer[i][j];
            }
        }
        memcpy(mxGetPr(mPWork), mPWorkMx, uWidth * bufSize * sizeof(real_T));
        mxSetFieldByNumber(simSnap, 0, 0, mPWork);
    }
    {
        mxArray* mIndWrt = mxCreateDoubleMatrix(1, uWidth, mxREAL);
        memcpy(mxGetPr(mIndWrt), idxWrt, uWidth * sizeof(int_T));
        mxSetFieldByNumber(simSnap, 0, 1, mIndWrt);
    }
    {
        mxArray* mBufSz = mxCreateDoubleMatrix(1, 1, mxREAL);
        memcpy(mxGetPr(mBufSz), &bufSize, sizeof(int_T));
        mxSetFieldByNumber(simSnap, 0, 2, mBufSz);
    }

    return simSnap;
}


/* Function: mdlSetSimState
*   Unpack the MATLAB structure passed and restore it into the RunTimeData
*   structure
*/
static void mdlSetSimState(SimStruct* S, const mxArray* simSnap) {
    int_T    i, j, bufSize;    
    real_T** buffer   = (real_T**)ssGetPWork(S)[PWORKS_IDX_BUFFER];
    int_T*   idxWrt   = (int_T*)ssGetPWork(S)[PWORKS_IDX_WRITE_IDX];
    int_T    uWidth  = (int_T)ssGetInputPortWidth(S, 0);
    {
        const mxArray* mBufSz = mxGetField(simSnap, 0, fieldNames[2]);
        memcpy(&bufSize, mxGetPr(mBufSz), sizeof(int_T));
        ssSetRWorkValue(S, IDX_BUFFER_SIZE, (real_T)bufSize);
    }
    {
        const mxArray* mPWork   = mxGetField(simSnap, 0, fieldNames[0]);
        real_T*        mPWorkMx = mxCalloc(uWidth * bufSize, sizeof(real_T));
        memcpy(mPWorkMx, mxGetPr(mPWork), uWidth * bufSize * sizeof(real_T));

        for (i = 0; i < uWidth; i++) {
            for (j = 0; j < bufSize; j++) {
                buffer[i][j] = mPWorkMx[i * bufSize + j];
            }
        }
    }
    {
        const mxArray* mIndWrt = mxGetField(simSnap, 0, fieldNames[1]);
        memcpy(idxWrt, mxGetPr(mIndWrt), uWidth * sizeof(int_T));
    }
}
#endif


#ifdef MATLAB_MEX_FILE
     #include "simulink.c"
#else
     #include <cg_sfun.h>
#endif
