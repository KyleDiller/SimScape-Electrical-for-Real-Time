function[MATE]=MATE_Sfun_DataFormat(rlc, source, switches, yout,Adp,B1dp,B2dp,Cdp,Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p,mateblk,u0,x0, statenames, rlcnames, sourcenames,switchtype,NumberOfNodes,PortType, Ts, SwitchVf)
% function d'analyse nodale pour determiner les parametres de calcul MATE


% statenames
% rlcnames
% sourcenames  To be added later

MATE.sizes=[size(Adp{1},1) size(B1dp{1},2) size(Ddp{1},1), size(switches,1) NumberOfNodes size(Adp,2) sum(PortType==1)]; %nb_state nb_input nb_output nb_switches nb_nodes nb_permutations nb_ItypePorts
MATE.rlc=rlc;
MATE.source=source;
MATE.yout=yout;
MATE.switches=switches;
MATE.Adp=Adp;
MATE.B1dp=B1dp;
MATE.B2dp=B2dp;
MATE.Cdp=Cdp;
MATE.Ddp=Ddp;
MATE.Yp=Yp;
MATE.Yindex=1:size(Yp{1},1);    %% assume the order is correct (1 2 3...) because we control the group
MATE.Cinj_p=Cinj_p;
MATE.Dinj_p=Dinj_p;
MATE.z_p=z_p;
MATE.x_p=x_p;
MATE.Ts=Ts;
MATE.nb_state=size(Adp{1},1);
MATE.nb_input=size(Ddp{1},2);
MATE.nb_output=size(Ddp{1},1);
MATE.nb_switch=size(switches,1);
MATE.x0=x0;
MATE.u0=u0;
MATE.NumberOfNodes=NumberOfNodes; % to be extended later
MATE.UinOrder=1:size(source,1);  % to be refined later SPS.OPALssn.colIN{i}
MATE.SwOrder=1:size(switches,1); %%%%%%%SPS.OPALssn.switchesIN{i}
%MATE.SwType= ones(size(switches,1),1); %% all ideal switches  SPS.OPALssn.switches_type{i}
MATE.SwType=switchtype;
%%% control of switch to be made externally in Simulink (ex: breaker)

%MATE.PortType=mateblk{1,3}; %to be refined later =MATEtype;
MATE.PortType=PortType;
% MATEtype=='V' or 'I'    %SPS.OPALssn.PortType_num{i};
if sum(MATE.PortType)==0
    MATE.PortType_num=0;
elseif sum(MATE.PortType)==length(MATE.PortType)
    MATE.PortType_num=1;
else
    MATE.PortType_num=9; % mixed type
end

MATE.Thy=[]; %%%%%SPS.OPALssn.FirstAndLastThy{i};
MATE.Thy_ind_src=[]; %%%SPS.OPALssn.Thy_ind_src(j)

MATE.rlcnames=rlcnames;
MATE.sourcenames=sourcenames;

% puting matrices into vectors for real-time

rowoffset=128;
MATE.rowoffset=rowoffset;

Adps_off=ceil(MATE.sizes(1)*MATE.sizes(1)/rowoffset)*rowoffset;  %each permutation starting memory location is rounded to the rowoffset
Adps=zeros(1,Adps_off*MATE.sizes(6));
for sw=1:MATE.sizes(6)
    for i=1:MATE.sizes(1)
        for j=1:MATE.sizes(1)
            Adps((sw-1)*Adps_off+(i-1)*MATE.sizes(1)+j)=Adp{sw}(i,j);
        end
    end
end
MATE.Adps=Adps;
MATE.SwitchVf=SwitchVf;






