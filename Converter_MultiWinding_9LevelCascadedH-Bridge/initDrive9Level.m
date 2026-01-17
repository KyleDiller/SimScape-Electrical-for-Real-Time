

% Source parameter
Src_V=3000;   % Voltage exprime  L-L RMS
Src_R=0.01;   % R source



% transfo parameters

T_Pbase= 500e3; % nominal power  Moteur 1MW  approx: transfo d'entree 1.5MVA
T_fbase=60;  % nominal frequency
T_Vpri=3000; % Primary voltage
T_Vpri_15comp=1/0.935; % facteur de compensation du transfo 15 deg
T_Vpri_15comp=1.1168;
T_Vsec=920;  % secondary voltages


Pbase=500e3;   % Transfo secondary impedance base
Zbase=T_Vsec*T_Vsec/T_Pbase;
wbase=2*pi*T_fbase;             % base de frequence angulaire
Lbase=Zbase/wbase;              % base d'inductance



Rtot_pu=0.01;
Ltot_pu=0.04;
T_2wind_ratio=0.2677;

%%%%  Lpu of secondary of 2 winding transfo
T_Rpu=Rtot_pu*9/10;  % 90% of R attributed to secondary in PU
T_Lpu=Ltot_pu*9/10;  % 90% of L attributed to secondary in PU
%%%%  Lpu of secondary of 3 winding transfo in series with 2e one (zig-zag)
T_Rpu2=Rtot_pu*5/10;  % 50% of R attributed to secondary in PU
T_Lpu2=Ltot_pu*5/10;  % 50% of L attributed to secondary in PU
%%%%  Lpu of primary
T_Lpu_pri=Ltot_pu/10;   %  10% attributed to primary
%disp('T_Lpu_primaire=0')
T_Rpu_pri=Rtot_pu/10;   

T_Lm_factor=inf;   % enlever inductance de magnetisation si =inf


    

Rgnd_rectifier=1e4; % ground referencing






%%%%%%%%%%%%%%%
% DC link capacitor
%%%%%%%%%%%%%

Cdclink=54e-3;
Rdclink_series=0;





