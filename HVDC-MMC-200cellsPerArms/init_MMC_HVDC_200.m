

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ts = 35e-6;	  % main time step
Tctrl=Ts;     % controller time step (must be=Ts) 


numberOfCell=200;
nb_cell=numberOfCell;
CarrierFreq=150; %Hz

 % Floating Capacitor value (F)
CC =0.005;
%IGBT Ron Roff
Ron = 0.01;%0.1
Roff=10000;


Ll = 0.04889; % MMC branch inductance value (H)
Rl= 0;          


E = 290e3;%100e3; % DC pole-ground voltage (V)
Rdc=3e6;



% AC System parameters

P_AC = 1000e6;%200e6; % Rated Active power (W)
V_AC_prim = 400e3;% % Rated Primary side rms L-L voltage (V)
V_AC_sec = 320e3; % Rated secondary side rms L-L voltage (V)

freq_left = 50; 
freq_right = 50; 

S_Xfo = P_AC*1.4; % Transformer MVA rating




