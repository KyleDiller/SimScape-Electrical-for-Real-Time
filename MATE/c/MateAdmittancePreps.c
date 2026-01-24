#include<stdio.h>
#include<stdlib.h>
#include "MatrixFun.h"
#include "MateUtils.h"
int MateAdmittancePreps(double Ad[], double Bd2[], double Cd[], double Dd[],double Yperm[],double Cinj_perm[], double Dinj_perm[], double z_perm[], double x_perm[], double B2inj_perm[], int NSTATE, int NINPUT, int NOUTPUT, int NNODE, int NITYPE)
{
    int not_alloc=-4;
    //nb_input_internal=size(Bd1,2)-nb_nodal_nodes;
    int nb_input_internal=NINPUT-NNODE;
    //Cinj_perm{i}=Cd(end-nb_nodal_nodes+1:end,:);
    //Dinj_perm{i}=Dd(end-nb_nodal_nodes+1:end,1:nb_input_internal);
    for(int i=0;i<NNODE;i++){
        for(int j=0;j<NSTATE;j++){
            Cinj_perm[i*NSTATE+j]=Cd[i*NSTATE+j+NSTATE*(NOUTPUT-NNODE)];
        }
        for(int j=0;j<nb_input_internal;j++){
            Dinj_perm[i*nb_input_internal+j]=Dd[i*NINPUT+j+NINPUT*(NOUTPUT-NNODE)];
        }
    }
    double *Cinj_permBU;
    double *Dinj_permBU;
    if ((Cinj_permBU=(double*)malloc((NNODE*NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((Dinj_permBU=(double*)malloc((NNODE*nb_input_internal+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    double *y;
    double *x1;
    double *x2;
    if ((y=(double*)malloc(((NNODE-NITYPE)*(NNODE-NITYPE)+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((x1=(double*)malloc(((NNODE-NITYPE)*(NITYPE)+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((x2=(double*)malloc(((NNODE-NITYPE)*(NITYPE)+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    double *Btmp;
    double *Ctmp;
    if ((Btmp=(double*)malloc((NNODE*NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    if ((Ctmp=(double*)malloc((NNODE*NSTATE+1)*sizeof(double)))  == NULL) {
		return(not_alloc);
    }
    double *tmp,*tmp1XX,*tmp2,*tmp3,*eye, *tmp4;
        if ((tmp=(double*)malloc((NITYPE*NNODE+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }
        if ((tmp1XX=(double*)malloc((NITYPE*(NNODE-NITYPE)+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }
        //printf("tmp1XX size %i\n", (NITYPE*(NNODE-NITYPE)+1) *sizeof(double));
        //printf("tmp1XX size %i addr: %p\n", (NITYPE*(NNODE-NITYPE)+1) *sizeof(double),(void *)tmp1XX);
        
        if ((tmp2=(double*)malloc((NITYPE*(NNODE-NITYPE)+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }
        if ((tmp3=(double*)malloc(((NNODE-NITYPE)*(NNODE-NITYPE)+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }
        if ((eye=(double*)malloc(((NNODE-NITYPE)*(NNODE-NITYPE)+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }
        if ((tmp4=(double*)malloc(((NNODE)*(NNODE)+1)*sizeof(double)))  == NULL) {
		    return(not_alloc);
        }

        // double *inv_z;
        // if ((inv_z=(double*)malloc(((NITYPE)*(NITYPE)+1)*sizeof(double)))  == NULL) {
		//     return(not_alloc);
        // }

    
    //z_perm{i}=[];
    //x_perm{i}=[];
    //%(C_no*B2_no+Dno_no)  Impedance ou admittance nodale discrete equivalente
    //Yperm{i}=Cd(end-nb_nodal_nodes+1:end,:)*Bd2(:,end-nb_nodal_nodes+1:end)+Dd(end-nb_nodal_nodes+1:end,end-nb_nodal_nodes+1:end) ;
    //          Yperm= Cd(X)
    //  <----NSTATE-------->
//      |---------|--------| _
//      |         |         NOUT-NNODE
//      |------------------| -
//      |   X     |  X     |
//      |         |        |NNODE
//      |---------|--------|
    //          Yperm= Bd2(X)
    //  <NIN-NNODE><-NNODE->
//      |---------|--------| _
//      |         |  X     |
//      |------------------| STATE
//      |         |  X     |
//      |         |        | -
//      |---------|--------|
    //          Yperm= +D(X)
    //  <NIN-NNODE><-NNODE->
//      |---------|--------| _
//      |         |         NOUT-NNODE
//      |------------------| -
//      |         |  X     |
//      |         |        |NNODE
//      |---------|--------|

        
    for(int i=0;i<NSTATE;i++){
        for(int j=0;j<NNODE;j++){
            Btmp[i*NNODE+j]=Bd2[i*NINPUT+j+(NINPUT-NNODE)];
        }
    }
    for(int i=0;i<NNODE;i++){
        for(int j=0;j<NSTATE;j++){
            Ctmp[i*NSTATE+j]=Cd[i*NSTATE+j+(NOUTPUT-NNODE)*NSTATE];
        }
    }
    mulmx(Yperm, Ctmp, Btmp, NNODE, NNODE, NSTATE);  // C*Bd2
    for(int i=0;i<NNODE;i++){
        for(int j=0;j<NNODE;j++){
            Yperm[i*NNODE+j]+=Dd[i*NINPUT+j+NINPUT*(NOUTPUT-NNODE)+(NINPUT-NNODE)];
        }  // C*Bd2+D
    }

    //printf("YpermTMP=");
    //printmat(Yperm, NNODE, NNODE); //OK
    int Ptype;
    if (NNODE==NITYPE) Ptype=1;
    else if (NITYPE==0) Ptype=0;
    else Ptype=2;
    
    if (Ptype==1){

        //Yperm{i}=inv(Yperm{i}); %4 oct 2017: revert to inv()
        matinv2(Yperm,NNODE);

        mulmx(Cinj_permBU, Yperm, Cinj_perm, NNODE, NSTATE, NNODE);
        mulmx(Dinj_permBU, Yperm, Dinj_perm, NNODE, nb_input_internal, NNODE);
//         void gain_mult_mx( double in1[], double gain,double out[], int row, int col)
// {
// 		for(int i=0;i<row*col;i++){
// 			out[i]=in1[i]*gain;  /* out=in1[]*in2 */
// 		}
// }

        gain_mult_mx(Cinj_permBU, -1.0 , Cinj_perm, NNODE, NSTATE);
        gain_mult_mx(Dinj_permBU, -1.0 , Dinj_perm, NNODE, nb_input_internal);
        // modif 23 sept.
        // cpmx(Cinj_perm, Cinj_permBU, NNODE*NSTATE);
        // cpmx(Dinj_perm, Dinj_permBU, NNODE*nb_input_internal);

        //Cinj_perm{i}=-Yperm{i}*Cinj_perm{i};
        //Dinj_perm{i}=-Yperm{i}*Dinj_perm{i};
        for (int i=0; i<NSTATE; i++){
            for (int j=0; j<NNODE; j++){
                B2inj_perm[i*NNODE+j]=Bd2[i*NINPUT+j +(NINPUT-NNODE)];
            }
        }
        
    }
    else if (Ptype==0){
        for (int i=0; i<NSTATE; i++){
            for (int j=0; j<NNODE; j++){
                B2inj_perm[i*NNODE+j]=Bd2[i*NINPUT+j +(NINPUT-NNODE)];
            }
        }
    }
    else{
        //j=1:(sum(portType==1));
        // j is set of 1:NITYPE
        //k=j(end)+1:length(portType);
        // so k is NITYPE+1:NNODE
        //z=Yperm{i}(j,j);
        double *inv_z;
        inv_z=z_perm;
        for(int i=0;i<NITYPE;i++){
            for(int j=0;j<NITYPE;j++){
                inv_z[i*NITYPE+j]=Yperm[i*NNODE+j];
            }
        }
        //% inv(z) is the only inverse in mixed V/I-type
        //% we replace inv(z) by inv_z computed by ldl() or chol()
        //%LDLTcode=max(LDLTcode,max(max(abs(z-z'))));   %test for symmetry
        //LDLTcode=min(LDLTcode,rcond(z));  %test for conditionning
        //inv_z=inv(z);   %4 oct 2017: revert to inv()
        matinv2(inv_z,NITYPE);
            //printf("inv_z=");
        //printmat(inv_z, NITYPE,NITYPE); //OK


//            Yperm
//      <----NNODE-------->
//      <-NITYPE->< ----  ><NNODE-NITYPE
//      |---------|--------| _
//      |   invz  |    x1  | NITYPE
//      |------------------| -
//      |         |        |
//      |   x2    |    y   | NNODE-NITYPE
//      |---------|--------| -
        //y=Yperm{i}(k,k);
        
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<(NNODE-NITYPE);j++){
                y[i*(NNODE-NITYPE)+j]=Yperm[NITYPE*NNODE+NITYPE+i*NNODE+j];
            }
        }
        //printf("y=");
        //printmat(y, (NNODE-NITYPE),(NNODE-NITYPE)); 
        
        //x1=Yperm{i}(j,k);
        for(int i=0;i<NITYPE;i++){
            for(int j=0;j<(NNODE-NITYPE);j++){
                x1[i*(NNODE-NITYPE)+j]=Yperm[NITYPE+i*NNODE+j];
            }
        }

        
        //x2=Yperm{i}(k,j);
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<NITYPE;j++){
                x2[i*NITYPE+j]=Yperm[NITYPE*NNODE+i*NNODE+j];
            }
        }

  
        //%YYY=Yperm{i}
        //%Yperm{i}=[ inv_z -inv_z*x1 ; x2*inv_z (y*z-x1*x2)*inv_z];
        //Yperm{i}=[ inv_z -inv_z*x1 ; x2*inv_z y-x2*inv_z*x1];
//      <----NNODE-------->
//      <-NITYPE->< ----  ><NNODE-NITYPE
//      |---------|--------| _
//      |   invz  |-invz*x1| NITYPE
//      |------------------| -
//      |         |        |
//      | x2*invz |        |<=y-x2*inv_z*x1
//      |---------|--------|
        //Yperm{i}=[ inv_z * ; * *];
        for(int i=0;i<NITYPE;i++){
            for(int j=0;j<NITYPE;j++){
                Yperm[i*NNODE+j]=inv_z[i*NITYPE+j];
            }
        }
        //Yperm{i}=[ * -inv_z*x1 ; * *];

 
        mulmx(tmp1XX, inv_z, x1, NITYPE, (NNODE-NITYPE), NITYPE);//tmp1XX=inv_z*x1  (I,I)*(I)(V)->(I,V)
        for(int i=0;i<NITYPE;i++){
            for(int j=0;j<(NNODE-NITYPE);j++){
                Yperm[NITYPE+i*NNODE+j]= tmp1XX[i*(NNODE-NITYPE)+j] *(-1.0);  // negate here
            }
        }
        
         //Yperm{i}=[ * * ; x2*inv_z *];
        mulmx(tmp2, x2, inv_z, (NNODE-NITYPE), NITYPE, NITYPE); //tmp2=x2*inv_z    (V,I)*(I,I)=(V,I)
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<NITYPE;j++){
                Yperm[NITYPE*NNODE+i*NNODE+j]=tmp2[i*(NITYPE)+j];  
            }
        }
        
        //Yperm{i}=[ * * ;  * y-x2*inv_z*x1];
        mulmx(tmp3, x2, tmp1XX, (NNODE-NITYPE),(NNODE-NITYPE), NITYPE); //x2*inv_z*x1    (V,I)*(I,V)->(V,V)
        submx(y, tmp3, y, (NNODE-NITYPE), (NNODE-NITYPE));             //y-x2*inv_z*x1    (V,V)-(V,V)=(V,V)
        //submx(y, y, tmp3, (NNODE-NITYPE), (NNODE-NITYPE));             //avoid crash
  
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<(NNODE-NITYPE);j++){
                Yperm[NITYPE*NNODE+NITYPE+i*NNODE+j]=y[i*(NNODE-NITYPE)+j];  
            }
        }
        
        //LDLTcode=max(LDLTcode,max(max(abs(Yperm{i}-Yperm{i}'))));  %test for symmetry




        //tmp=[-inv_z;-x2*inv_z];
        //tmp1XX=[zeros(length(k),length(j)) eye(length(k))];
        //tmp=[tmp tmp1XX'];
        // Cinj_perm NNODE  NSTATE
        // Dinj_perm NNODE  nb_input_internal
//      <----NNODE-------->
//      <-NITYPE->< ---- - ><NNODE-NITYPE
//      |---------|--------| _
//      |   -invz |    0   | NITYPE
//      |------------------| -
//      |         |        |--
//      |-x2*inv_z|    1   | NNODE-NITYPE
//      |---------|--------|--
        if (1==0){
        for(int i=0;i<NNODE;i++){
            for(int j=0;j<NNODE;j++){
                tmp4[i*NNODE+j]=0.0;  // initialize to 0
            }
        }

        //for(int i=0;i<NNODE;i++){
        for(int i=0;i<NITYPE;i++){
            for(int j=0;j<NITYPE;j++){
                tmp4[i*NNODE+j]=0.0-Yperm[i*NNODE+j];  // negate of NITYPE cols of Yperm
            }
        }
        ident(eye, (NNODE-NITYPE));  // identify matrix
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<(NNODE-NITYPE);j++){
                tmp4[(NNODE*NITYPE)+(NITYPE)+i*NNODE+j]=eye[i*(NNODE-NITYPE)+j];
            }
        }

        gain_mult_mx(tmp2, -1.0 , tmp2, NITYPE ,(NNODE-NITYPE));
        for(int i=0;i<(NNODE-NITYPE);i++){
            for(int j=0;j<(NITYPE);j++){
                tmp4[(NNODE*NITYPE)+i*NNODE+j]=tmp2[i*(NITYPE)+j];
            }
        }
        }
        else {
            for(int i=0;i<NITYPE;i++){
                for(int j=0;j<(NNODE-NITYPE);j++){
                    tmp4[NITYPE+i*NNODE+j]=0.0;  // initialize to 0
                }
            }
            for(int i=0;i<NNODE;i++){
                for(int j=0;j<NITYPE;j++){
                    tmp4[i*NNODE+j]=0.0-Yperm[i*NNODE+j];  //
                }
            }
            ident(eye, (NNODE-NITYPE));  // identify matrix
            for(int i=0;i<(NNODE-NITYPE);i++){
                for(int j=0;j<(NNODE-NITYPE);j++){
                    tmp4[(NNODE*NITYPE)+(NITYPE)+i*NNODE+j]=eye[i*(NNODE-NITYPE)+j];
                }
            }

        }


        //WHERE IS  -x2*inv_z  !!!!???????

        //%Cinj_perm{i}=[-inv(z) 0*eye(length(k));-x2*inv(z) eye(length(k))]*Cinj_perm{i};
        //%Dinj_perm{i}=[-inv(z) 0*eye(length(k));-x2*inv(z) eye(length(k))]*Dinj_perm{i};
        //Cinj_perm{i}=tmp*Cinj_perm{i};
        //Dinj_perm{i}=tmp*Dinj_perm{i};
    
        mulmx(Cinj_permBU, tmp4, Cinj_perm, NNODE, NSTATE, NNODE);

        mulmx(Dinj_permBU, tmp4, Dinj_perm, NNODE, (NINPUT-NNODE), NNODE);

        cpmx(Cinj_perm, Cinj_permBU, NNODE*NSTATE);
        cpmx(Dinj_perm, Dinj_permBU, NNODE*(NINPUT-NNODE));

        
        //z_perm{i}=-inv_z;
        gain_mult_mx( inv_z, -1.0 , inv_z, NITYPE,NITYPE);
        cpmx(z_perm, inv_z, NITYPE*NITYPE);

        //x_perm{i}=x1;
        cpmx(x_perm, x1, (NNODE-NITYPE)*NITYPE);
        
        for (int i=0; i<NSTATE; i++){
            for (int j=0; j<NNODE; j++){
                B2inj_perm[i*NNODE+j]=Bd2[i*NINPUT+j +(NINPUT-NNODE)];
            }
        }

    }
        free(y);
        free(x1);
        free(x2);
        //printf("tmp1XX size %i addr: %p\n", (NITYPE*(NNODE-NITYPE)+1) *sizeof(double),(void *)tmp1XX);

        free(tmp1XX);
        free(tmp);
        free(tmp2);
        free(tmp3); 
        free(eye);
        free(Btmp);
        free(Ctmp);
        free(tmp4);
        free(Cinj_permBU);
        free(Dinj_permBU);
        //free(inv_z);
}
    
        
    
    
