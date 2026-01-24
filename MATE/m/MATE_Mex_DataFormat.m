function[MATERT]=MATE_Mex_DataFormat(Ac,Bc,Cc,Dc, switches,SwType,disc,u0,x0, Ts, sizes,SwitchVf)

As=[];
Bs=[];
for i=1:size(Ac,1)      % put in series to avoid col-row issues with matlab fortran C...
    As=[As Ac(i,:)];
    Bs=[Bs Bc(i,:)];
end
Cs=[];
Ds=[];
for i=1:size(Cc,1)
    Cs=[Cs Cc(i,:)];
    Ds=[Ds Dc(i,:)];
end

MATERT.As= As;  
MATERT.Bs= Bs;
MATERT.Cs= Cs;
MATERT.Ds= Ds;
if size(switches,1)>0
    MATERT.Rsw=switches(:,4)';
else
    MATERT.Rsw=[];
end
MATERT.SwType=SwType;
MATERT.h= Ts;
MATERT.DISC= disc;
MATERT.U0= u0;
MATERT.X0= x0;
MATERT.sizes=sizes;
MATERT.SwitchVf=SwitchVf;

% MATERT.sizes,MATERT.As,MATERT.Bs,MATERT.Cs,MATERT.Ds,MATERT.Rsw,MATERT.SwType,MATERT.h, MATERT.DISC,MATERT.U0,MATERT.X0





