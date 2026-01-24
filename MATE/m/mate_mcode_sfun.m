function mate_mcode_sfun(block,x)
%MSFUNTMPL_BASIC A template for a Leve-2 M-file S-function
%   The M-file S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the M-file S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2007 The MathWorks, Inc.

%%
%% The setup method is used to setup the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction


%% Function: setup ===================================================
%% Abstract:
%%   Set up the S-function block's basic characteristics such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

MATE=block.DialogPrm(1).Data;
% Register number of ports

block.NumInputPorts  = 3;  % Simulink source input for ABCD equation, Nodal voltages, switch gates
    
block.NumOutputPorts = 4;  % Y, Ihist, all ABCD output, switch status




% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

Nw_input=size(MATE.Ddp{1},2); % -MATE.NumberOfNodes;
Nb_output=size(MATE.Ddp{1},1);



% Override input port properties
% V-I source input
block.InputPort(1).Dimensions        = max([ 1 Nw_input]);
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'Real';
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(1).SamplingMode      = 0;

% nodal voltage
block.InputPort(2).Dimensions        = MATE.NumberOfNodes;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'Real';
block.InputPort(2).DirectFeedthrough = true;
block.InputPort(2).SamplingMode      = 0;

% switch gates
block.InputPort(3).Dimensions        = max([size(MATE.switches,1) 1]);
block.InputPort(3).DatatypeID        = 0;  % double
block.InputPort(3).Complexity        = 'Real';
block.InputPort(3).DirectFeedthrough = true;
block.InputPort(3).SamplingMode      = 0;


% Override output port properties Y
block.OutputPort(1).Dimensions       = [MATE.NumberOfNodes,MATE.NumberOfNodes];
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'Real';
block.OutputPort(1).SamplingMode     = 0;

% Override output port properties Ihist
block.OutputPort(2).Dimensions       = [MATE.NumberOfNodes];
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 0;

% Override output port properties Ihist
block.OutputPort(2).Dimensions       = [MATE.NumberOfNodes];
block.OutputPort(2).DatatypeID       = 0; % double
block.OutputPort(2).Complexity       = 'Real';
block.OutputPort(2).SamplingMode     = 0;


block.OutputPort(3).Dimensions       = Nb_output;
block.OutputPort(3).DatatypeID       = 0; % double
block.OutputPort(3).Complexity       = 'Real';
block.OutputPort(3).SamplingMode     = 0;

block.OutputPort(4).Dimensions       = max([size(MATE.switches,1) 1]);
block.OutputPort(4).DatatypeID       = 0; % double
block.OutputPort(4).Complexity       = 'Real';
block.OutputPort(4).SamplingMode     = 0;



% Register parameters
block.NumDialogPrms     = 1;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [MATE.Ts 0];    % take sample time of group 1

%% -----------------------------------------------------------------
%% The M-file S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
% % % % SPS=block.DialogPrm(1).Data;
MATE=block.DialogPrm(1).Data;



block.NumDworks = 11;   %

% les variables d'etats ou de status sont stockées ci-apres
% % % % for i=1:ssn_groups
i=1;
  % % % % % nb_state=SPS.OPALssn.dims{i}(1);
  % % % % % nb_input=SPS.OPALssn.dims{i}(2);
  % % % % % nb_output=SPS.OPALssn.dims{i}(3);
nb_state=MATE.nb_state;
nb_input=MATE.nb_input;
nb_output=MATE.nb_output;
nb_switch=MATE.nb_switch;
nb_nodal_nodes=MATE.NumberOfNodes;


  block.Dwork(1).Name            = ['xn_g' num2str(i)];    % states of partitions
  block.Dwork(1).Dimensions      = nb_state + isequal(nb_state,0);
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;
  
  block.Dwork(2).Name            = ['un_g' num2str(i)];   % inputs of partitions  
  block.Dwork(2).Dimensions      = nb_input;
  block.Dwork(2).DatatypeID      = 0;      % double
  block.Dwork(2).Complexity      = 'Real'; % real
  block.Dwork(2).UsedAsDiscState = true;
  
  block.Dwork(3).Name            = ['sw_g' num2str(i)];   % switch status of partitions  
  block.Dwork(3).Dimensions      = nb_switch+ isequal(nb_switch,0);
  block.Dwork(3).DatatypeID      = 0;      % double
  block.Dwork(3).Complexity      = 'Real'; % real
  block.Dwork(3).UsedAsDiscState = true;
  
  block.Dwork(4).Name            = ['sw_volta' num2str(i)];   % switch voltage of partitions  
  block.Dwork(4).Dimensions      = nb_switch+ isequal(nb_switch,0);
  block.Dwork(4).DatatypeID      = 0;      % double
  block.Dwork(4).Complexity      = 'Real'; % real
  block.Dwork(4).UsedAsDiscState = true;
  
  block.Dwork(5).Name            = ['sw_volta_old' num2str(i)];   % switch voltage (previous step) of partitions  
  block.Dwork(5).Dimensions      = nb_switch+ isequal(nb_switch,0);
  block.Dwork(5).DatatypeID      = 0;      % double
  block.Dwork(5).Complexity      = 'Real'; % real
  block.Dwork(5).UsedAsDiscState = true;
  
  block.Dwork(6).Name            = ['sw_sta_old_g' num2str(i)];   % switch status of partitions   
  block.Dwork(6).Dimensions      = nb_switch+ isequal(nb_switch,0);
  block.Dwork(6).DatatypeID      = 0;      % double
  block.Dwork(6).Complexity      = 'Real'; % real
  block.Dwork(6).UsedAsDiscState = true;
  
  block.Dwork(7).Name            = ['Ihist' num2str(i)];   % I_history   
  block.Dwork(7).Dimensions      = nb_nodal_nodes;
  block.Dwork(7).DatatypeID      = 0;      % double
  block.Dwork(7).Complexity      = 'Real'; % real
  block.Dwork(7).UsedAsDiscState = true;
  
  block.Dwork(8).Name            = ['Yold' num2str(i)];   % admittance Y  
  block.Dwork(8).Dimensions      = nb_nodal_nodes*nb_nodal_nodes;
  block.Dwork(8).DatatypeID      = 0;      % double
  block.Dwork(8).Complexity      = 'Real'; % real
  block.Dwork(8).UsedAsDiscState = true;
  
  block.Dwork(9).Name            = ['scrtmp' num2str(i)];   % old_source of partitions  
  block.Dwork(9).Dimensions      = nb_input;
  block.Dwork(9).DatatypeID      = 0;      % double
  block.Dwork(9).Complexity      = 'Real'; % real
  block.Dwork(9).UsedAsDiscState = true;
  
  block.Dwork(10).Name            = ['xtmp' num2str(i)];   %states 
  block.Dwork(10).Dimensions      = nb_state + isequal(nb_state,0);
  block.Dwork(10).DatatypeID      = 0;      % double
  block.Dwork(10).Complexity      = 'Real'; % real
  block.Dwork(10).UsedAsDiscState = true;

  block.Dwork(11).Name            = ['vf' num2str(i)];   %switch offset 
  block.Dwork(11).Dimensions      = nb_switch + isequal(nb_switch,0);
  block.Dwork(11).DatatypeID      = 0;      % double
  block.Dwork(11).Complexity      = 'Real'; % real
  block.Dwork(11).UsedAsDiscState = true;



%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)
MATE=block.DialogPrm(1).Data;


    if isempty(MATE.switches)
        x=[0];
        y=-1;
    else
        x=MATE.switches(:,3)';
        y=x*0-1;
    end
    if ~isempty(MATE.Adp{1})
        block.Dwork(1).Data=MATE.x0;
    else
        block.Dwork(1).Data=0;
    end
    block.Dwork(2).Data=0*MATE.u0;
    block.Dwork(7).Data=zeros(1,MATE.NumberOfNodes);
    block.Dwork(8).Data=zeros(1,MATE.NumberOfNodes*MATE.NumberOfNodes);
    
     block.Dwork(3).Data=x*0;
     block.Dwork(4).Data=y+1;
     block.Dwork(5).Data=y*0;
     if ~isempty(MATE.Adp{1})
        block.Dwork(10).Data=0*MATE.x0; 
     end
     %vf=block.DialogPrm(2).Data;
     %vf=MATE.SwitchVf;
     if MATE.nb_switch>0
        block.Dwork(11).Data=MATE.SwitchVf; 
     else
         block.Dwork(11).Data=0;
     end
    % src_old=block.Dwork(2).Data; 
    % I_histold=block.Dwork(7).Data;
    % Yvec=block.Dwork(8).Data;

%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this 
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
function Start(block)

%block.Dwork(1).Data(1) = 0;

%endfunction



%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%

function Outputs(block)
MATE=block.DialogPrm(1).Data;
a=[1 2 3];
%block.OutputPort(1).Data = block.Dwork(1).Data + block.InputPort(1).Data+block.Dwork(2).Data+sum(a) ;
%         [  x          ]n   [u_sw..  u_in   .. u_no]n     [u_sw..  u_in   .. u_no]n+1
%                                                          [        u             ]n+1
%  
%[ ]      [             ]    [                      ]      [             |        ]      
%[x]n+1 = [     A       ] +  [         B1           ]   +  [    B2_in    |  B2_no ]
%[ ]      [             ]    [                      ]      [             |        ]
%
%         [  x          ]n+1
%[    ]   [             ]                                  [             |        ]
%[yin ]   [    Cin      ]                                  [    Din_in   | Din_no ]  <- Din
%[    ]   [             ]                                  [             |        ]
%[ -- ]=  [ ---------   ]                               +  [  ------------------- ] -----------
%[    ]   [             ]                                  [             |        ]
%[yno ]   [    Cno      ]                                  [    Dno_in   | Dno_no ]
%[ n+1]   [             ]                                  [             |        ]

% x[n+1]=A*x[n]+B1*u[n]+B2_in*u_in[n+1]   +  B2_no*u_no[n+1]
% y_no[n+1]=C_no*(A*x[n]+B1*u[n]+B2_in*u_in[n+1]   +  B2_no*u_no[n+1]) +Dno_in*u_in[n+1] +Dno_no*u_no[n+1]
% y_no[n+1]=C_no*(A*x[n]+B1*u[n]+B2_in*u_in[n+1])+ Dno_in*u_in[n+1]   + (C_no*B2_no+Dno_no)*u_no[n+1])
%                  Term from History                                              nodal Z or Y 
% y_in[n+1]=Cin*x[n+1]+Din*u[n+1]  

%% step 1: read nodal voltage and complete the state-equations


Uin=block.InputPort(1).Data;   % source signals coming from Simulink (commanded I or V source from SPS)

V=block.InputPort(2).Data;   % Nodal voltage 

SWin=block.InputPort(3).Data;   % switch gate signals from Simulink 

MATE=block.DialogPrm(1).Data;

block.NumDworks = 11;   %  number of vector to be saved

% completion of group equations with the V solution from nodal problem YV=I.

  nb_state=MATE.nb_state;
  nb_input=MATE.nb_input;
  nb_nodal_nodes=MATE.NumberOfNodes;
  nb_input_internal=nb_input-nb_nodal_nodes;
  nb_output=MATE.nb_output;
  nb_switch=MATE.nb_switch;
  Yindex=MATE.Yindex;

  first_last_thyristor=MATE.Thy; % to be refined later

  % get last switch status from the last time step.
  sw=0;
  for j=1:nb_switch
      if block.Dwork(3).Data(j)>0.5
          sw=sw+2^(nb_switch-j);
      end
  end
sw=sw+1;   % 0 logic to 1 logic


    A=MATE.Adp{sw};
    B=MATE.B1dp{sw};
    B2=MATE.B2dp{sw};
    C=MATE.Cdp{sw};
    D=MATE.Ddp{sw};
    Y=MATE.Yp{sw};
    Cinj=MATE.Cinj_p{sw};
    Dinj=MATE.Dinj_p{sw};
    zz=MATE.z_p{sw};   
    xx=MATE.x_p{sw};    

    x_tmp_old=block.Dwork(10).Data; 
    src_old=block.Dwork(2).Data; 
    I_histold=block.Dwork(7).Data;
    Yvec=block.Dwork(8).Data;
    for kkk=1:nb_nodal_nodes
       Yold(kkk,:)=Yvec((kkk-1)*nb_nodal_nodes+1:(kkk-1)*nb_nodal_nodes+nb_nodal_nodes);
    end
   % Yold=Y; %test
    % block.Dwork(10).Data=x_tmp; % save for later completion of state calcujlation at the start of the next time step
    % block.Dwork(7).Data=I_hist;
    % block.Dwork(2).Data=UIn_group_ordered;

    %%% should not need Yindex=SPS.OPALssn.NodeIndex{i};
    PortType=MATE.PortType; 
    NumPortTypeVolt=sum(PortType==0);

        if PortType==0
            % completion of state calculation with newly available V
            if nb_state>0.5
                block.Dwork(1).Data=x_tmp_old+B2(:,end-nb_nodal_nodes+1:end)*V(Yindex);
            end

            % concatenation of input [internal nodal]
            Vall=[src_old(1:nb_input_internal); V(Yindex)];

            % disp([ 'vallm=' num2str(Vall')])
            % disp(['xm=' num2str(block.Dwork(1).Data')])


        elseif PortType==1

            tmp=Yold*V(Yindex)+I_histold;  % the multiplication of Ihist by Y is now made in pre-calculation (23oct 2009) Minus sign made in fts6perm (3 Nov)
            Vall=[src_old(1:nb_input_internal); tmp];
            if nb_state>0.5
                block.Dwork(1).Data=x_tmp_old+B2(:,nb_input_internal+1:end)*tmp;
            end
             % disp([ 'MMMM start'])
             % disp([ 'vallm=' num2str(Vall')])
             % disp([ 'temp1=' num2str(tmp')])
             % disp([ 'x_tmp_old=' num2str(x_tmp_old')])
             % xxxx=block.Dwork(1).Data(1:6);
             % disp(['xm=' num2str(xxxx')]);
             % disp(['Ihistold=' num2str(I_histold')]);

        else  %mixed type
            % Voltage input is easy. Voltage are ordered last.

            V1=V(Yindex);
            V2=V1(end-NumPortTypeVolt+1:end);
            V3=V1(1:end-NumPortTypeVolt);
             %      Ityp Vtyp
             %  Vx=[ tmp   V2 ]
            tmp=zz*(xx*V2-V3)+I_histold(1:end-NumPortTypeVolt);
            if nb_state>0.5
                block.Dwork(1).Data=x_tmp_old+B2(:,end-nb_nodal_nodes+1:end)*[tmp;V2];
            end
            % concatenation of input [internal nodal]
            Vall=[src_old(1:nb_input_internal); tmp ; V2];
        end

        % sw=sw
        % TEMP1=(xx*V2-V3)'
        % XX=xx
        % V2=V2
        % V3=V3


    
    % calculation of  outputs
    if nb_state>0.5
        Yout=C*block.Dwork(1).Data+D*Vall;
    else
        Yout=D*Vall;
    end
    
    %block.Dwork(8).Data=Yout;
    block.OutputPort(3).Data=Yout;  % output Yout;


    if nb_switch>0.5
        block.Dwork(5).Data=block.Dwork(4).Data;
        block.Dwork(4).Data=Yout(1:nb_switch);
    end



    if nb_switch>0.5
       block.OutputPort(4).Data=block.Dwork(3).Data;  % output switch status
    end


 %% step 2 :compute the next switch status, next Y, next Ihist
  
  % copy old switch status
  block.Dwork(6).Data=block.Dwork(3).Data;

  UIn_group_ordered=Uin(MATE.UinOrder);   % source 

  % save the source input for source vector completion at the start of the
  % next time step
  block.Dwork(9).Data=UIn_group_ordered;

  if nb_switch>0.5
     SWInOrder=SWin(MATE.SwOrder);   % switch signals 
  end
 
  SwType=MATE.SwType;
  
 % selection de la permutation ABCD du aux switchs



 sw=0;

 for j=1:nb_switch

     switch SwType(j)
         case 1  %switch
             if SWInOrder(j)>0.5
                 block.Dwork(3).Data(j)=1;
             else
                 block.Dwork(3).Data(j)=0;
             end            
         case 2  % breaker
             if SWInOrder(j)>0.5
                 block.Dwork(3).Data(j)=1;
             else if ~isequal(sign(block.Dwork(4).Data(j)),sign(block.Dwork(5).Data(j)))
                       block.Dwork(3).Data(j)=0;
                  end
             end
         case {3 , 32} %diode
             pred=2*block.Dwork(4).Data(j)-block.Dwork(5).Data(j);
             xxx=pred;
             % line change: with iteration method it's better not to make
             % prediction, it change cause switch toggling at each time
             % step see case SSN_iterative_power_hvdc_itvc2.mdl
             pred=block.Dwork(4).Data(j);

             %if (block.Dwork(3).Data(j)>0 & block.Dwork(2).Data(j)==0) | (pred>0 & (block.Dwork(2).Data(j)==1)) 
             if block.Dwork(4).Data(j)<0
                    block.Dwork(3).Data(j)=1;
             else
                    block.Dwork(3).Data(j)=0;
             end

         case 4  %thyristor
             pred2=2*block.Dwork(4).Data(j)-block.Dwork(5).Data(j);
             pred=block.Dwork(4).Data(j);
             %if (block.Dwork(3).Data(j)>0)  & ( (SWInOrder(j)>0.000001)| (block.Dwork(2).Data(j)==1))
             %if ((block.Dwork(4).Data(j)>0)  & (SWInOrder(j)>0.000001)) | (pred>0 & (block.Dwork(3).Data(j)==1))
             if (SWInOrder(j)>0.5 & j==2)
                 xxxxx=1;
             end
             %%%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%
             %%%% Vsw<-Vf !!!!
             %%%%%%%%%%%%%%%%%
             if pred<-block.Dwork(11).Data(j)  & ( (SWInOrder(j)>0.5) | (block.Dwork(3).Data(j)>0.5))
             %if (pred>0 & block.Dwork(3).Data(j)<0.5 & SWInOrder(j)>0.5) | (pred>0 & (block.Dwork(3).Data(j)>0.5)) 
                 % % % if (block.Dwork(3).Data(j)==0)
                 % % %     % compensation IVC ici
                 % % %     %if SPS.OPALssn.Params.ivc==1
                 % % %     if (SWInOrder(j))<(1-1e-9)
                 % % %         %SWInOrder(j)
                 % % %        %UIn_group_ordered{i}(SPS.OPALssn.Thy_ind_src{i}(j))=-block.Dwork(3).Data(j)*(SWInOrder(j)-1); 
                 % % %        UIn_group_ordered(MATE.Thy_ind_src)=-block.Dwork(4).Data(j)*(SWInOrder(j)); 
                 % % %        %qwer=SWInOrder(j)  % test for fts5conversion bug
                 % % %     else
                 % % %     end
                 % % % end
                 block.Dwork(3).Data(j)=1;         
             else
                 block.Dwork(3).Data(j)=0;
                 
             end
         case {5 , 6}    % GTO MOSFET IGBT 
             %disp([num2str(block.Dwork(3).Data(j))  num2str(SWInOrder(j))])
             if (block.Dwork(4).Data(j)>0)  &  (SWInOrder(j)>0.00001)
                 block.Dwork(3).Data(j)=1;         
             else
                 block.Dwork(3).Data(j)=0;
             end           
     end

 end



for j=1:nb_switch
      if block.Dwork(3).Data(j)>0.5
         %sw=sw+2^(j-1);
         sw=sw+2^(nb_switch-j);
      end
end

sw=sw+1;   % 0 logic to 1 logic

if sw==1
    qqqq=1;
end
if sw==8
    qqqq=2;
end


 % end SWITCH SECTION
 

    block.Dwork(6).Data(1)=0; 

 
  PortType=MATE.PortType;  

  
  % x[n+1]=A*x[n]+B1*u[n]+B2_in*u_in[n+1]

          
    A=MATE.Adp{sw};
    B=MATE.B1dp{sw};
    B2=MATE.B2dp{sw};
    C=MATE.Cdp{sw};
    D=MATE.Ddp{sw};
    Y=MATE.Yp{sw};
    Cinj=MATE.Cinj_p{sw};
    Dinj=MATE.Dinj_p{sw};
    
    %UIn_group_ordered(nb_switch+1:nb_input_internal,:)=UIn_group_ordered(nb_switch+1:nb_input_internal,:).* block.Dwork(3).Data;
    UIn_group_ordered(1:nb_switch,:)=UIn_group_ordered(1:nb_switch,:).* block.Dwork(3).Data;
  
    if nb_state>0.5

        x_tmp=A*block.Dwork(1).Data ...;
            +B*Vall ...;  % this is all source including nodal solution known from past history
            +B2(:,1:nb_input_internal)*UIn_group_ordered(1:nb_input_internal,:);  %  this is only internal source are known at the current step.

        I_hist=Cinj*x_tmp+Dinj*UIn_group_ordered(1:nb_input_internal,:);
    else
        x_tmp=[];
        I_hist= Dinj*UIn_group_ordered(1:nb_input_internal,:);
    end

    % XTMP=x_tmp'
    % IHIST=I_hist'
    % % VALL=Vall'
    %          disp([ 'xtmp' num2str(x_tmp')])
    %          xxxx=block.Dwork(1).Data(1:6);
    %          disp(['x=' num2str(xxxx')]);
    %          disp(['Ihistold(end=' num2str(I_hist')]);
    %          disp([ 'MMMM end'])


  % y_no[n+1]=C_no*(A*x[n]+B1*u[n]+B2_in*u_in[n+1])+ Dno_in*u_in[n+1] 

 
% compute Y
bigY=zeros(nb_nodal_nodes,nb_nodal_nodes); % reset to 0
bigI=zeros(nb_nodal_nodes,1); % reset to 0
Yvec=[];
for kkk=1:nb_nodal_nodes
     Yvec=[Yvec Y(kkk,:)];
end

% 
% for kkk=1:length(Yindex)
%     bigI(Yindex(kkk))=bigI(Yindex(kkk))-I_hist(kkk);
% end


 
block.OutputPort(1).Data=Y;  % output Y
block.OutputPort(2).Data=-I_hist;  % output Ihist

block.Dwork(8).Data=Yvec;
if nb_state>0
    block.Dwork(10).Data=x_tmp; % save for later completion of state calcujlation at the start of the next time step
end
block.Dwork(7).Data=I_hist;
block.Dwork(2).Data=UIn_group_ordered;



  
%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
function Update(block)

%block.Dwork(1).Data(1) = block.InputPort(1).Data;

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%end Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)


%end Terminate




