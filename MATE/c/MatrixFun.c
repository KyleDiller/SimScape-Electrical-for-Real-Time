#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define SQRT2 1.414213562373095
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

int putmx(double *A,double *W, int a_c, int w_r, int w_c, int a_off){

	for(int i=0;i<w_r;i++){
		for(int j=0;j<w_c;j++){
			A[a_off+i*a_c+j]=*W++;
		}
	}
	return(0);
}


void cpmx(double out[], double in[], int n)
{

	for(int i=0;i<n;i++){
		out[i]=in[i];
	}

}
void transposemx(double out[], double in[], int nrow, int ncol)
{
    double *outptr=out;
	for(int i=0;i<nrow;i++){
		for(int j=0;j<ncol;j++){
			out[j*nrow+i]=in[j+i*ncol];
		}
	}

}




void mulmx(double out[], double in1[], double in2[], int row, int col, int irc)
{

	double accu;

   	for(int i=0;i<row;i++){
		for(int j=0;j<col;j++){
			accu=0.0;
			for (int k=0;k<irc;k++){
				accu=accu+in1[i*irc+k]*in2[j+k*col];
			}
			out[i*col+j]=accu;
		}
	}
}

void ident(double I[], int n){

	for (int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			I[i*n+j]=0.0;
			if (i==j){
				I[i*n+j]=1.0;
			}
		}
	}
}
void opmx(double out[], double in1[], double in2[], int row, int col,int opt)
{


	if (opt==0){  /* addition */
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]+in2[i];  /* out=in1+in2 */
		}
	}


	if (opt==1){  /* soustraction */
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]-in2[i];    /* out=in1-in2 */

		}
	}

	if (opt==2){  /* multiplication matrice*scalaire */
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]*in2[0];  /* out=in1[]*in2 */
		}
	}
}

void addmx(double in1[], double in2[], double out[], int row, int col)
{
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]+in2[i];  /* in1+in2 =>out*/
		}
}


void submx(double in1[], double in2[], double out[], int row, int col)
{
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]-in2[i];    /* in1-in2=>out */

		}
}

void gain_mult_mx( double in1[], double gain,double out[], int row, int col)
{
		for(int i=0;i<row*col;i++){
			out[i]=in1[i]*gain;  /* out=in1[]*in2 */
		}
}


/* Matrix inversion ------------------------------------------- */
int matinv2(double *A,int N)
{
    int            icol,irow,i,j,k,l,ll;
    double         big,dum,pivinv;
	/* int            ipiv[10000],indxr[10000],indxc[10000];   */
	int *ipiv,*indxr,*indxc;
    int            singular=0;
	int	not_alloc=-9;
    icol = 0;
    irow = 0;
	if ((ipiv=(int*)malloc(N*sizeof(int)))  == NULL) {
		return(not_alloc);
    }
	if ((indxr=(int*)malloc(N*sizeof(int))) == NULL) {
		return(not_alloc);
    }
	if ((indxc=(int*)malloc(N*sizeof(int))) == NULL) {
		return(not_alloc);
    }

	for (j=0;j<N;j++) ipiv[j]=0;
    for (i=0;i<N;i++) {
		big=0.0;
		for (j=0;j<N;j++)
			if (ipiv[j] != 1)
				for (k=0;k<N;k++) {
					if (ipiv[k] == 0) {
						if (fabs(A[j*N+k]) >= big) {
							big=fabs(A[j*N+k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) singular=1;

				}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l=0;l<N;l++){
						SWAP(A[irow*N+l],A[icol*N+l])
					}

				}
				indxr[i]=irow;
				indxc[i]=icol;
				if (A[icol*N+icol] == 0.0) singular=2;

				pivinv=1.0/A[icol*N+icol];
				A[icol*N+icol]=1.0;
				for (l=0;l<N;l++) A[icol*N+l] *= pivinv;

				for (ll=0;ll<N;ll++)
					if (ll != icol) {
						dum=A[ll*N+icol];
						A[ll*N+icol]=0.0;
						for (l=0;l<N;l++) A[ll*N+l] -= A[icol*N+l]*dum;
					}
    }
    for (l=N-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<N;k++)
                SWAP(A[k*N+indxr[l]],A[k*N+indxc[l]]);
    }
	free(ipiv);
	free(indxr);
	free(indxc);
    return singular;
}

void mulmxcpx(double out_r[], double out_i[], double in1_r[], double in1_i[],int conj1, double in2_r[], double in2_i[], int conj2, int row, int col, int irc)
{
	int i,j,k;    /* out=in1*in2 */
	double accu_r,accu_i;
	double sign1=1.0;
	double sign2=1.0;
	if (conj1==1) sign1=-1.0;
	if (conj2==1) sign2=-1.0;

   	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			accu_r=0.0;
			accu_i=0.0;
			for (k=0;k<irc;k++){
				accu_r=accu_r+in1_r[i*irc+k]*in2_r[j+k*col]-sign1*in1_i[i*irc+k]*sign2*in2_i[j+k*col];
				accu_i=accu_i+in1_r[i*irc+k]*sign2*in2_i[j+k*col]+sign1*in1_i[i*irc+k]*in2_r[j+k*col];
			}
			out_r[i*col+j]=accu_r;
			out_i[i*col+j]=accu_i;
		}
	}
}

/* COMPLEX Matrix inversion ------------------------------------------- */
int matinvcpx(double *Ar, double *Ai, int N)
{
    int            icol,irow,i,j,k,l,ll;
    double         big;

	int *ipiv,*indxr,*indxc;
	double pivinv_r,pivinv_i,*Artmp,dum_r,dum_i;
    int            singular=0;
	int	not_alloc=-9;
    icol = 0;
    irow = 0;

	if ((ipiv=(int*)malloc(N*sizeof(int)))  == NULL) {
		return(not_alloc);
    }
	if ((indxr=(int*)malloc(N*sizeof(int))) == NULL) {
		return(not_alloc);
    }
	if ((indxc=(int*)malloc(N*sizeof(int))) == NULL) {
		return(not_alloc);
    }
	if ((Artmp=(double*)malloc(N*N*sizeof(double))) == NULL) {
		return(not_alloc);
    }
	for (j=0;j<N;j++) ipiv[j]=0;
    for (i=0;i<N;i++) {
		big=0.0;
		for (j=0;j<N;j++)
			if (ipiv[j] != 1)
				for (k=0;k<N;k++) {
					if (ipiv[k] == 0) {
						// if (fabs(A[j*N+k]) >= big) {
						if (Ar[j*N+k]*Ar[j*N+k]+Ai[j*N+k]*Ai[j*N+k] >= big) {
							//big=fabs(A[j*N+k]);
							big=Ar[j*N+k]*Ar[j*N+k]+Ai[j*N+k]*Ai[j*N+k];
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) singular=1;

				}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l=0;l<N;l++){
						//SWAP(A[irow*N+l],A[icol*N+l])
						SWAP(Ar[irow*N+l],Ar[icol*N+l])
							SWAP(Ai[irow*N+l],Ai[icol*N+l])
					}

				}
				indxr[i]=irow;
				indxc[i]=icol;
				//if (A[icol*N+icol] == 0.0) singular=2;
				if ((fabs(Ar[icol*N+icol])+ fabs(Ai[icol*N+icol]))== 0.0) singular=2;
				pivinv_r=Ar[icol*N+icol]/(Ar[icol*N+icol]*Ar[icol*N+icol]+Ai[icol*N+icol]*Ai[icol*N+icol]);
				pivinv_i=(0.0-Ai[icol*N+icol])/(Ar[icol*N+icol]*Ar[icol*N+icol]+Ai[icol*N+icol]*Ai[icol*N+icol]);
				//pivinv=1.0/A[icol*N+icol];
				//A[icol*N+icol]=1.0;
				Ar[icol*N+icol]=1.0;
				Ai[icol*N+icol]=0.0;

				//for (l=0;l<N;l++) A[icol*N+l] *= pivinv;
				for (l=0;l<N;l++) Artmp[icol*N+l] = Ar[icol*N+l]*pivinv_r - Ai[icol*N+l]*pivinv_i;
				for (l=0;l<N;l++) Ai[icol*N+l] = Ai[icol*N+l]*pivinv_r + Ar[icol*N+l]*pivinv_i;
				for (l=0;l<N;l++) Ar[icol*N+l] = Artmp[icol*N+l];

				for (ll=0;ll<N;ll++)
					if (ll != icol) {
						//dum=A[ll*N+icol];
						dum_r=Ar[ll*N+icol];
						dum_i=Ai[ll*N+icol];
						//A[ll*N+icol]=0.0;
						Ar[ll*N+icol]=0.0;
						Ai[ll*N+icol]=0.0;
						//for (l=0;l<N;l++) A[ll*N+l] -= A[icol*N+l]*dum;
						for (l=0;l<N;l++) Artmp[ll*N+l] = Ar[ll*N+l] - Ar[icol*N+l]*dum_r + Ai[icol*N+l]*dum_i;
						for (l=0;l<N;l++) Ai[ll*N+l] = Ai[ll*N+l] - Ar[icol*N+l]*dum_i - Ai[icol*N+l]*dum_r;
						for (l=0;l<N;l++) Ar[ll*N+l] = Artmp[ll*N+l];
					}
    }
    for (l=N-1;l>=0;l--) {
        if (indxr[l] != indxc[l])
            for (k=0;k<N;k++) {
                //SWAP(A[k*N+indxr[l]],A[k*N+indxc[l]]);
                SWAP(Ar[k*N+indxr[l]],Ar[k*N+indxc[l]]);
                SWAP(Ai[k*N+indxr[l]],Ai[k*N+indxc[l]]);
            }
    }
	free(ipiv);
	free(indxr);
	free(indxc);
	free(Artmp);
    return singular;
}

