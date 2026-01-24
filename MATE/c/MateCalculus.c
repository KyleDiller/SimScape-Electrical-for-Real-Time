#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "MatrixFun.h"
#define MAX_NB_SWITCH 20
int MateDiscretize(double A[], double B1[],double B2[], double C[], double D[], double h, int opt, int holdtype, int NSTA, int NIN, int NOUT)
{

	double a1,a2,a3,a4;
	double *I;
	int not_alloc=-9;
	double *tmp,*tmp1,*denx,*Ad,*Bd1,*Bd2,*Ah,*Ctmp,*Dtmp,*Bf=NULL,*Bz=NULL;
	double fa=(2.0-sqrt(2.0))/2.0;
	double *tmp0_r, *tmp0_i, *tmp1_r, *tmp1_i, *den1x_r, *den1x_i, *den2x, *num1_r, *num1_i;

	int jz,jf,j,i,k;
	denx = 0;
	tmp = 0;
	Ad=0;
	Bd1=0;


	if ((I=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
	ident(I,NSTA);
	if ((Ah=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }


	


	switch (opt) {

	case 5:			/* Radau IIA(5) */
		if ((den1x_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((den1x_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((den2x=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((num1_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((num1_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp0_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp0_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1_r=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1_i=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		//den_gain=-0.016666666666666666666666666666666666666666666666667;
		//den1= (A*h - (2.681082873627752133895790743 + 3.0504301992474105694263776247*i)*I);
		//den1t=(A*h - (2.681082873627752133895790743 - 3.0504301992474105694263776247*i)*I);
		//den2= (A*h - 3.6378342527444957322084185135777*I)    ;

		//num_gain=0.05000;
		//num1=(A*h + (4. + 2.00*i)*I);
		//num1t=(A*h + (4. - 2.00*i)*I);
		//Ad=real(num_gain/den_gain*inv(den1)*num1t*inv(den1t)*num1*inv(den2))

        a1=-2.681082873627752133895790743;
        a2=-3.0504301992474105694263776247;
        a3=- 3.6378342527444957322084185135777;

        opmx(Ah,A,&h,NSTA,NSTA,2);          //A*h

        opmx(tmp0_r,I,&a1,NSTA,NSTA,2);        //-2.681082873627752133895790743*I
        opmx(den1x_r,Ah,tmp0_r,NSTA,NSTA,0);   // A*h -2.681082873627752133895790743*I
        opmx(den1x_i,I,&a2,NSTA,NSTA,2);    // -3.0504301992474105694263776247*i*I
        matinvcpx(den1x_r,den1x_i,NSTA);    //inv(A*h -2.681082873627752133895790743*I -3.0504301992474105694263776247*i*I)

        opmx(tmp0_r,I,&a3,NSTA,NSTA,2);        // -3.6378342527444957322084185135777*I
        opmx(den2x,Ah,tmp0_r,NSTA,NSTA,0);     //(A*h - 3.6378342527444957322084185135777*I)
        matinv2(den2x,NSTA);                // inv(A*h - 3.6378342527444957322084185135777*I)

        a1= 4.0;
        a2= 2.0;
        a3=-3.0;

        opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
        opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h +4*I

        opmx(num1_i,I,&a2,NSTA,NSTA,2);     // 2*i*I

        mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,1,NSTA,NSTA,NSTA);
        mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
        mulmxcpx(tmp0_r,tmp0_i,tmp1_r,tmp1_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
        mulmx(tmp1_r,tmp0_r,den2x,NSTA,NSTA,NSTA);
        opmx(A,tmp1_r,&a3,NSTA,NSTA,2);


		if (holdtype==1){
			// 1/2*I - 2/15*h*A + 1/60*h*h*A*A=
			//.01666666666666666666666666666666666666667 *
			// (A*h - (4. + 3.741657386773941385583748732316549301756*i)*I)*
			// (A*h - (4. - 3.741657386773941385583748732316549301756*i)*I)


			a1=-4.0;
			a2=-3.741657386773941385583748732316549301756;
			a3=-1.0*h;

			opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
			opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h - 4*I

			opmx(num1_i,I,&a2,NSTA,NSTA,2);     // -3.741657386773941385583748732316549301756*i*I

			mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
			mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
			mulmxcpx(tmp0_r,tmp0_i,tmp1_r,tmp1_i,0,num1_r,num1_i,1,NSTA,NSTA,NSTA);
			mulmx(tmp1_r,tmp0_r,den2x,NSTA,NSTA,NSTA);
			opmx(tmp1_r,tmp1_r,&a3,NSTA,NSTA,2);


			mulmx(Bd2,tmp1_r,B1,NSTA,NIN,NSTA);	/* Bd2 not Bd1!   */

			// 1/2*I + 1/30*h*A =
			// 0.03333333333333*(A*h + 15*I)

			a1=15.0;
			a2=-2.0*h;

			opmx(tmp0_r,I,&a1,NSTA,NSTA,2);
			opmx(num1_r,Ah,tmp0_r,NSTA,NSTA,0);    // A*h +15*I
			opmx(num1_i,I,I,NSTA,NSTA,1);    //  0*i*I

			mulmxcpx(tmp0_r,tmp0_i,den1x_r,den1x_i,0,num1_r,num1_i,0,NSTA,NSTA,NSTA);
			mulmxcpx(tmp1_r,tmp1_i,tmp0_r,tmp0_i,0,den1x_r,den1x_i,1,NSTA,NSTA,NSTA);
			mulmx(tmp0_r,tmp1_r,den2x,NSTA,NSTA,NSTA);
			opmx(tmp1_r,tmp0_r,&a2,NSTA,NSTA,2);

			mulmx(Bd1,tmp1_r,B1,NSTA,NIN,NSTA);	/* Bd1 , not Bd2  ! */



			
			cpmx(B1,Bd1,NSTA*NIN);
			cpmx(B2,Bd2,NSTA*NIN);
		
		}
		
		free(tmp0_r);
		free(tmp0_i);
		free(tmp1_r);
		free(tmp1_i);
		free(den1x_r);
		free(den1x_i);
		free(den2x);
		free(num1_r);
		free(num1_i);
		free(Bd1);
		free(Bd2);
		break;


	case 31:
		fa=(2.0+sqrt(2.0))/2.0;
	case 3:
		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp1=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a4=0.0-fa;
		opmx(Ah,A,&h,NSTA,NSTA,2);          // A*h
		opmx(denx,Ah,&a4,NSTA,NSTA,2);     // -fa*h*A
		opmx(denx,I,denx,NSTA,NSTA,0);	// I-fa*h*A
		matinv2(denx,NSTA);


		a1=1.0-2.0*fa;
		opmx(tmp,Ah,&a1,NSTA,NSTA,2);       //  (1-2*fa)*h*A
		opmx(tmp,I,tmp,NSTA,NSTA,0);	    // I + (1-2*fa)*h*A

		mulmx(tmp1,denx,tmp,NSTA,NSTA,NSTA);// (I+(1-2*fa)*h*A)*denx;
		mulmx(A,tmp1,denx,NSTA,NSTA,NSTA);// Ad=(I+(1-2*fa)*h*A)*denx*denx;

		if (holdtype==1){
			a1=2.0*fa-fa*fa;
			a2=0.0-fa*fa;
			opmx(tmp,Ah,&a2,NSTA,NSTA,2);  // -fa*fa*h*A
			opmx(tmp1,I,&a1,NSTA,NSTA,2);  //   (2*fa-fa*fa)*I
			opmx(tmp,tmp,tmp1,NSTA,NSTA,0);	//  (2*fa-fa*fa)*I- fa*fa*h*A

			mulmx(tmp1,denx,tmp,NSTA,NSTA,NSTA);
			mulmx(tmp,tmp1,denx,NSTA,NSTA,NSTA);
			opmx(tmp1,tmp,&h,NSTA,NSTA,2);
			//mulmx(Bd1,tmp1,B1,NSTA,NIN,NSTA);	/*  Bd1=((2*fa-fa*fa)*I-fa*fa*h*A)*denx*h*B  */

		    mulmx(Bd1,tmp1,B1,NSTA,NIN,NSTA);	/* Bd1   */

			a3=(1.0-2.0*fa+fa*fa)*h;
			mulmx(tmp,denx,denx,NSTA,NSTA,NSTA);
			opmx(tmp,tmp,&a3,NSTA,NSTA,2);

			//mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);	/*  	Bd2=(1-2*fa+fa*fa)*denx*h*B;  */

			mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);	/* Bd2   */

			
			
			cpmx(B1,Bd1,NSTA*NIN);
			cpmx(B2,Bd2,NSTA*NIN);
			
			//cpmx(B1,Bd1,NSTA*NIN);
			//cpmx(B2,Bd2,NSTA*NIN);
		}
		

		free(denx);
		free(tmp);
		free(tmp1);
		free(Bd1);
		free(Bd2);

		break;

	case 2:

		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a1=0.0-0.5;
		a2=0.5;
		opmx(Ah,A,&h,NSTA,NSTA,2);  /* A*h */
		opmx(denx,Ah,&a1,NSTA,NSTA,2); /* -0.5*h*A */
		opmx(denx,I,denx,NSTA,NSTA,0);	/* I-0.5*h*A */
		matinv2(denx,NSTA);
		opmx(tmp,Ah,&a2,NSTA,NSTA,2); /* 0.5*h*A */
		opmx(tmp,I,tmp,NSTA,NSTA,0);	/* I+0.5*h*A */
		mulmx(A,denx,tmp,NSTA,NSTA,NSTA);  /* 	Ad=denx*(I+h*A/2) */


		opmx(tmp,denx,&h,NSTA,NSTA,2); /* denx*h */
		mulmx(Bd1,tmp,B1,NSTA,NIN,NSTA);  /* denx*h*B	 */
		opmx(Bd2,Bd1,&a2,NSTA,NIN,2); /* denx*h*B/2 */
		cpmx(B1,Bd2,NSTA*NIN);
	    cpmx(B2,Bd2,NSTA*NIN); 
	    
		free(denx);
		free(tmp);
		free(Bd1);
		free(Bd2);


		break;
		
	case 1:

		if ((denx=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
        if ((tmp=(double*)malloc((NSTA*NSTA+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		if ((Bd1=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
		}
		if ((Bd2=(double*)malloc((NSTA*NIN+1)*sizeof(double)))  == NULL) {
			return(not_alloc);
        }
		a1=1.0;
		a2=0.0;
		opmx(Ah,A,&h,NSTA,NSTA,2);  /* A*h */
		opmx(denx,I,Ah,NSTA,NSTA,1);	/* I-h*A */
		matinv2(denx,NSTA); // denx=inv(I-hA)
		cpmx(A,denx,NSTA*NSTA);
	
		opmx(tmp,A,&h,NSTA,NSTA,2); /* denx*h */
		mulmx(Bd2,tmp,B1,NSTA,NIN,NSTA);  /* denx*h*B	 */
		opmx(Bd1,Bd2,&a2,NSTA,NIN,2); /* Bd1=0 */
		cpmx(B1,Bd1,NSTA*NIN);
	    cpmx(B2,Bd2,NSTA*NIN); 
	    
		free(denx);
		free(tmp);
		free(Bd1);
		free(Bd2);


		break;
		

	}



	free(Ah);
	free(I);
	return(0);
}

int MatePermutize(double Amm[], double Bmm[], double Cmm[], double Dmm[], double Rswitch[],int permutation, int NSTATE, int NINPUT, int NOUTPUT, int NSWITCH)
{
// function [Adperm, Bdperm, Cdperm, Ddperm, Yperm, Cinj_perm, Dinj_perm, z_perm, x_perm, Bd2perm]= MATESwitchPermute(Am,Bm,Cm,Dm,Rswitch,Ts,nb_nodal_nodes,portType,disc)
// % FTS2PERM 	Calculate continuous time combinations of state-space matrices due to the switches

// IA=eye(length(Am));
double *IA, *I, Yswitch[MAX_NB_SWITCH],YswitchDiagonal[MAX_NB_SWITCH];
int ETATSW[MAX_NB_SWITCH], NIV_SWITCH[MAX_NB_SWITCH];
int not_alloc=-9;

	

if ((IA=(double*)malloc((NSTATE*NSTATE+1)*sizeof(double)))  == NULL) {
	return(not_alloc);
}
ident(IA,NSTATE);



// [NSTATE,NINPUT]=size(Bm); ok
// [NOUTPUT,NSTATE]=size(Cm); ok
// NSWITCH=length(Rswitch); ok
// ETATSW=-ones(NSWITCH,1);

for (int i=0; i<NSWITCH;i++) ETATSW[i]=-1;
    
// NIV_SWITCH=1:NSWITCH;

for (int i=0; i<NSWITCH;i++) NIV_SWITCH[i]=i;
// D2=zeros(NINPUT,NOUTPUT); not used
//	if ((D2=(int*)malloc((NINPUT*NOUTPUT+1)*sizeof(int)))  == NULL) {
//		return(not_alloc);
//    }

// Yswitch=1./Rswitch;

for (int i=0; i<NSWITCH;i++) Yswitch[i]=1.0/Rswitch[i];
    
// I=diag(ones(1,NOUTPUT));
if ((I=(double*)malloc((NOUTPUT*NOUTPUT+1)*sizeof(double)))  == NULL) {
	return(not_alloc);
}
ident(I,NOUTPUT);


// for i=1:2^NSWITCH

int idx;

int i;

i=permutation;
    for (int j=0; j<NSWITCH; j++) ETATSW[j]=0;
    int tmp=0;
    int div=1;
    for (int j=0; j<NSWITCH-1; j++) div=div*2;
    idx=i;
    for (int j=0; j<NSWITCH; j++) {
        if ( (idx-div)>=0) {
            ETATSW[j]=1;
            idx=idx-div;
        }
            div=div/2;
    }

    //printf("permuation=%i e3 %i  e2  %i  e1 %i \n", i,ETATSW[2],ETATSW[1],ETATSW[0]);



//     disp(['permutation ' num2str(i) ' of ' num2str(2^NSWITCH)]);
//     xstr=dec2bin(i-1,NSWITCH);
//     for j=1:NSWITCH
//         ETATSW(j)=isequal(xstr(j),'1');
//     end

//     %ETATSW=ones(1,NSWITCH);
//     if i>1
//         % sens des permutations	sw3: derniere switch de la liste
//         % perm	sw1 sw2	sw3
//         % 1: 	off off	off
//         % 2:	    off	off	on
//         % 3:  	off	on 	off
//         %	...
//         % 8	    on  on	on


//         % The following code is inspired from official SPS R14 code.
//         % This is acknoledged here.
//         % This is not a copy from Opal-RT code but inspired from SPS TLC code


//         N_SwitchOn=find(ETATSW==1);
//         YswitchDiagonal=zeros(NSWITCH,1);
//         YswitchDiagonal(N_SwitchOn)=Yswitch(N_SwitchOn);
    for (int i=0; i<NSWITCH; i++){ 
        if (ETATSW[i]==1) YswitchDiagonal[i]=Yswitch[i];
        else YswitchDiagonal[i]=0;
    }

    //cpmx(double out[], double in[], int n)

//         Amm=Am;  % travailler sur des copies de matrices d'origines car il y a une boucle sur les switch
//         Bmm=Bm;
//         Cmm=Cm;
//         Dmm=Dm;
    double temp, *DxCol, *BdCol, *tmp1, *tmp2;
    if ((DxCol=(double*)malloc((NOUTPUT+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((BdCol=(double*)malloc((NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((tmp1=(double*)malloc((NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((tmp2=(double*)malloc((NINPUT+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
  
    for(int k=0;k<NSWITCH;k++){

//         for k=1:length(ETATSW)
//             %temp = 1/(1-D[nSw*(nInputs+1)]*yswitch[nSw]*kSw);
//             temp=1/(1-Dmm(k,k)*YswitchDiagonal(k));
        temp=1.0/(1.0-Dmm[k*NINPUT+k]*YswitchDiagonal[k]);
//             %for(i=0; i<nOutputs; i++)  DxCol[i]=D[i*nInputs+nSw]*yswitch[nSw]*temp*kSw;
//             DxCol=Dmm(:,k)*YswitchDiagonal(k)*temp;
        for(int m=0;m<NOUTPUT;m++) DxCol[m]=Dmm[m*NINPUT+k]*YswitchDiagonal[k]*temp;
//             %DxCol[nSw] = temp;
//             DxCol(k,1)=temp;
        DxCol[k]=temp;
 



//             %for(i=0; i<nStates; i++) BDcol[i]=B[i*nInputs + nSw]*yswitch[nSw]*kSw;
//             BdCol=Bmm(:,k)*YswitchDiagonal(k);
        for(int m=0;m<NSTATE;m++) BdCol[m]=Bmm[m*NINPUT+k]*YswitchDiagonal[k];

//             %/* Copy row nSw of C into tmp1 and zero it out in C */
//             tmp1=Cmm(k,:);
//             Cmm(k,:)=tmp1*0;
        for(int m=0;m<NSTATE;m++) {
            tmp1[m]=Cmm[k*NSTATE+m];
            Cmm[k*NSTATE+m]=0;
        }
 

//             %  /* Copy row nSw of D into tmp2 and zero it out in D */
//             tmp2=Dmm(k,:);
//             Dmm(k,:)=tmp2*0;
        for(int m=0;m<NINPUT;m++) {
            tmp2[m]=Dmm[k*NINPUT+m];
            Dmm[k*NINPUT+m]=0;
        }


//             % C = C + DxCol * tmp1, D = D + DxCol * tmp2
//             Cmm= Cmm +DxCol*tmp1;
//             Dmm= Dmm +DxCol*tmp2;
        for(int m=0;m<NOUTPUT;m++) {
            for(int n=0;n<NSTATE;n++) Cmm[m*NSTATE+n]=Cmm[m*NSTATE+n]+DxCol[m]*tmp1[n];
            for(int n=0;n<NINPUT;n++) Dmm[m*NINPUT+n]=Dmm[m*NINPUT+n]+DxCol[m]*tmp2[n];
        }
 

//             % A = A + BdCol*C(nSw,:), B = B + BdCol*D(nSw,:)

//             Amm= Amm + BdCol* Cmm(k,:);
//             Bmm= Bmm + BdCol* Dmm(k,:);
        for(int m=0;m<NSTATE;m++) {
            for(int n=0;n<NSTATE;n++) Amm[m*NSTATE+n]=Amm[m*NSTATE+n]+BdCol[m]*Cmm[k*NSTATE+n];
            for(int n=0;n<NINPUT;n++) Bmm[m*NINPUT+n]=Bmm[m*NINPUT+n]+BdCol[m]*Dmm[k*NINPUT+n];
        }

//         end

//         if 1==1
//             As=Amm;
//             Bs=Bmm;
//             Cs=Cmm;
//             Ds=Dmm;
//         end

//     else
//         As=Am;
//         Bs=Bm;
//         Cs=Cm;
//         Ds=Dm;
//     end


    }

    free(tmp1);
    free(tmp2);
    free(BdCol);
    free(DxCol);
    
    free(I);
    free(IA);
    return(0);
}