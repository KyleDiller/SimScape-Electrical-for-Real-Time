int putmx(double *A,double *W, int a_c, int w_r, int w_c, int a_off);

void cpmx(double out[], double in[], int n);

void mulmx(double out[], double in1[], double in2[], int row, int col, int irc);

void ident(double I[], int n);

void opmx(double out[], double in1[], double in2[], int row, int col,int opt);

void transposemx(double out[], double in[], int nrow, int ncol);

void addmx(double in1[], double in2[], double out[], int row, int col);

void submx(double in1[], double in2[], double out[], int row, int col);

void gain_mult_mx( double in1[], double gain ,double out[], int row, int col);

/* Matrix inversion ------------------------------------------- */
int matinv2(double *A,int N);

void mulmxcpx(double out_r[], double out_i[], double in1_r[], double in1_i[],int conj1, double in2_r[], double in2_i[], int conj2, int row, int col, int irc);


/* COMPLEX Matrix inversion ------------------------------------------- */
int matinvcpx(double *Ar, double *Ai, int N);


