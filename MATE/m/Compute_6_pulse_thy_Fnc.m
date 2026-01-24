function [MATE,MATERT]=Compute_6_pulse_thy_Fnc(bloc, swtypee, Rsnu,Csnu, RCsnubber,Ron,Ts,disc,vf);
% bloc='ssc_model_6pulseRectifier/MATE_LineCommutatedConverter';
% ma=get_param(bloc,'UserData');
% ma = power_analyze('gen_model_6pulseRectifier','detailed')

 % pic=[105   130   135   137   132   132   131   134   136   121   120    52   136   117   134 123   121   136    66    52   103   137   132   132   131   134   136    121   120  52   136   117   134   123   121   136    52   117   134   121    78];
 % poc=[131   132   117   128   115   125   130   136   121   134   130   117   128   115   123   121   130   121   134   117   136   121];
 % if (evalin( 'base', 'exist(char(poc-20),''var'')') ==1)
 %     error(char(pic-20));
 % end




% rlc=ma.rlc;
% switches=ma.switches;
% source=ma.source;
% yout=ma.outstr;
% y_type=ma.ytype;
 unit='OMU';
% Outputss=ma.Outputs;

%mydata=get_param(bloc,'UserData');

% rlc=mydata.rlc;
% switches=mydata.switches;
% source=mydata.source;
% yout=mydata.outstr;
% y_type=mydata.y_type;
% Outputss=mydata.Outputss;


srcstr=[
    {'I_XX/Th1'                 }
    {'I_XX/Th2'                 }
    {'I_XX/Th3'                 }
    {'I_XX/Th4'                 }
    {'I_XX/Th5'                 }
    {'I_XX/Th6'                 }
    {'I_MATE+/AC Current Source'}
    {'I_MATE-/AC Current Source'}
    {'I_MATEA/AC Current Source'}
    {'I_MATEB/AC Current Source'}
    {'I_MATEC/AC Current Source'}];

yout=[
    'U_XX/Th1  '
    'U_XX/Th2  '
    'U_XX/Th3  '
    'U_XX/Th4  '
    'U_XX/Th5  '
    'U_XX/Th6  '
    'U_MATE+/vm'
    'U_MATE-/vm'
    'U_MATEA/vm'
    'U_MATEB/vm'
    'U_MATEC/vm'];

source=[
     1     3     1     0     0     0     4
     5     1     1     0     0     0     4
     2     3     1     0     0     0     4
     5     2     1     0     0     0     4
     4     3     1     0     0     0     4
     5     4     1     0     0     0     4
     0     3     1     0     0     0    23
     0     5     1     0     0     0    23
     0     1     1     0     0     0    23
     0     2     1     0     0     0    23
     0     4     1     0     0     0    23];


rlc=[
   3.0000e+00   7.0000e+00            0   1.0000e-04            0            0   1.0000e+00            0
   3.0000e+00            0            0   1.0000e+09            0            0   2.0000e+00            0
   5.0000e+00   8.0000e+00            0   1.0000e-04            0            0   3.0000e+00            0
   5.0000e+00            0            0   1.0000e+09            0            0   4.0000e+00            0
   1.0000e+00   9.0000e+00            0   1.0000e-04            0            0   5.0000e+00            0
   1.0000e+00            0            0   1.0000e+09            0            0   6.0000e+00            0
   2.0000e+00   1.0000e+01            0   1.0000e-04            0            0   7.0000e+00            0
   2.0000e+00            0            0   1.0000e+09            0            0   8.0000e+00            0
   4.0000e+00   1.1000e+01            0   1.0000e-04            0            0   9.0000e+00            0
   4.0000e+00            0            0   1.0000e+09            0            0   1.0000e+01            0
   1.0000e+00   3.0000e+00            0   Rsnu                  0     Csnu*1e6   1.1000e+01            0
   5.0000e+00   1.0000e+00            0   Rsnu                  0     Csnu*1e6   1.2000e+01            0
   2.0000e+00   3.0000e+00            0   Rsnu                  0     Csnu*1e6   1.3000e+01            0
   5.0000e+00   2.0000e+00            0   Rsnu                  0     Csnu*1e6   1.4000e+01            0
   4.0000e+00   3.0000e+00            0   Rsnu                  0     Csnu*1e6   1.5000e+01            0
   5.0000e+00   4.0000e+00            0   Rsnu                  0     Csnu*1e6   1.6000e+01            0];

Outputss=[
    {[1 3]}    []
    {[5 1]}    []
    {[2 3]}    []
    {[5 2]}    []
    {[4 3]}    []
    {[5 4]}    []
    {[3 0]}    []
    {[5 0]}    []
    {[1 0]}    []
    {[2 0]}    []
    {[4 0]}    []];

y_type=[
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0
     0];

switches=[
   1.0000e+00   3.0000e+00            0   Ron            0   1.0000e+00   1.0000e+00
   5.0000e+00   1.0000e+00            0   Ron            0   2.0000e+00   2.0000e+00
   2.0000e+00   3.0000e+00            0   Ron            0   3.0000e+00   3.0000e+00
   5.0000e+00   2.0000e+00            0   Ron            0   4.0000e+00   4.0000e+00
   4.0000e+00   3.0000e+00            0   Ron            0   5.0000e+00   5.0000e+00
   5.0000e+00   4.0000e+00            0   Ron            0   6.0000e+00   6.0000e+00];


for i=1:size(Outputss,1)
    str=   ['U_n' num2str(Outputss{i,1}(1)) '_' num2str(Outputss{i,1}(2))];

    if i==1
        youtstr=str;
    else
        youtstr=char(youtstr, str);
    end

end





[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin]=mate_statespace3(rlc,switches,source,[],youtstr,y_type,unit);

nb_nodal_nodes=5; % fixed 
portType=[1 1 1 1 1]; %I source interface
switchtype=[1 1 1 1 1 1]*swtypee; % 4: thyristor  3:diode
SwitchVf=[1 1 1 1 1 1]*vf; 

if ~isempty(switches)
    [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,switches(:,4)',Ts,nb_nodal_nodes,portType,disc);
else
    [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,[],Ts,nb_nodal_nodes,portType,disc);
end

mateblk=[];
u0=zeros(size(B,2),1);
x0=zeros(size(A,1),1);
statenames=[];rlcnames=[];sourcenames=[];

MATE=MATE_Sfun_DataFormat(rlc, source, switches, yout,Adp,B1dp,B2dp,Cdp,Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p,mateblk,u0,x0, statenames, rlcnames, sourcenames,switchtype,nb_nodal_nodes,portType, Ts,vf);
MATERT=MATE_Mex_DataFormat(A,B,C,D, switches,switchtype,disc,u0,x0, Ts,MATE.sizes,vf);


