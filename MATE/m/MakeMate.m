function [MATE, MATERT]= MakeMate(modele, disc)

%Ts=25e-6;
SILENT=0;  % 0 = less disp()


% 23 juillet  works on RTScape_netlist23juillet.slx
%zip('backup.zip', {'*.m', '*.mat'});
if 1==0
modele= 'RouteTestModel_1';   % multiple ground connected
modele= 'RouteTestModel_2';   % multiple ground not connected 
modele= 'RouteTestModel_3';   % multiple ground not connected 
modele= 'RouteTestModel_3A';  % model that make sense, with breakers
modele= 'RouteTestModel_4';  % model that make sense, with breakers, with m-file MATE Sfun
modele= 'RouteTestModel_4_a_Prep';  % simple current source ok
modele= 'RouteTestModel_4_b_Prep';  % simple voltage RL source
modele= 'RouteTestModel_4Prep';  % good!!!! after 4_a 4_b debug!!!
modele= 'RouteTestModel_6Prep';  %xfo 1ph built as Y-Y  : seems ok!!!
modele= 'RouteTestModel_7Prep_a';  %pi-line input (Norton case) ok
modele= 'RouteTestModel_7Prep';  %pi-line input +switch+xfo (Norton case) OKK!!
modele= 'RouteTestModel_7Prep_b';  %check Y (yout can have wrong sign!) OK
modele= 'RouteTestModel_8Prep_b';  %pi-line input +xfo (mix Norton=Thevenin case, 2 MATE interfaces) OK with unlimited custom variable resistance!!!
% with trap OK but NOT with art5 or art5-be!!!
% added mixed saturation to avoid +- eps for R
% also to note that the connection to the various interface is not so
% clear.... So the order is really [I-source V-source]  is case of mixed Norton-Thevenin 
% replace  'RouteTestModel_8Prep';   because include breakers.
modele= 'RouteTestModel_8Prep_e';  %Mix-works art5
modele= 'RouteTestModel_8Prep_d';  %mix works art5  %--->> we need to accept negative diagonla values!!!
modele= 'RouteTestModel_8Prep_b';  %19 aout: now works with art5: must allow disgonla negative values!!!

%modele= 'RouteTestModel_8Prep_c';   % with only Thevenin MATE +xfo+switch (correction dans la polarite des sources internes p/r a node1 node2!!!
%modele= 'RouteTestModel_9Prep';   % multiple Thevenin MATE. xfo breaker, internal source  GOOD 17 aout
%modele= 'RouteTestModel_10Prep';   % testing mix-MATE 1 ph  _b Vtype _c
%Itype   OK!
modele= 'RouteTestModel_11Prep'; modele_run='RouteTestModel_11';
modele= 'RouteTestModel_12_a/MATEPREP'; % 9 nodes mix
modele= 'RouteTestModel_13/MATEPREP';   % 6 node mix (simpler)  works in C!!!
modele= 'Route_MMC/MATEPREP'; 

end


model_top=get_param(modele,'Parent');

if isempty(model_top)
    Ts=evalin('base', get_param(bdroot,'FixedStep'));
else
    Ts=evalin('base', get_param(model_top,'FixedStep'));
end

%% DISCRETISE
if 1==0
disc=1; %1 backeuler  2 trap  5 pade(2,3)   3: pade(2,3_+BE_input
end

%modele=bdroot;
disp(['Running node analysis on ' modele]);

% at the start comment both MATE source
mateblk=[find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE Node')];
matetmp={};
for i=1:size(mateblk,1)
    try 
        a=get_param(mateblk{i},'conty');
        % if does not have this property, this is a single pphase block
    catch
        x=size(matetmp,1);
        matetmp{x+1,1}= mateblk{i};
    end
end
mateblk_sld=matetmp;   % a patch until we modify all test case with SLD MATE blocks with conty parameter

mateblk_1ph=[find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','MaskType','MATE Node','conty', 'Single-Phase')];



MakeMATEInterface([],mateblk_sld,mateblk_1ph,[],[],[],[],[],[],'commentMATE');

blk=find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','ReferenceBlock','ee_lib');
blk=[blk;find_system(modele, 'regexp', 'on', 'FollowLinks', 'on','LookUnderMasks','all','ReferenceBlock','fl_lib')]; % base ssc lib blocks


% easy way to find MATE type from the first block encountered...
%[mateblk] = MakeMATEInterface(blk,mateblk,[],[],[],'step1');






nb_blk=size(blk,1);

han=zeros(nb_blk,3);
rlc=zeros(1,6);
src=zeros(1,6);
sensors=zeros(1,6);
blkmask={
'Electrical Reference',[1],[1,0],-1; 
 ['Grounded Neutral' 10 '(Three-Phase)'],[1 ],[1,0],-1;
'Resistor',[1 2],[1 1],-1;                              % block type, electric ports indice for LConn and RConn in this order
'Inductor',[1 2],[1 1],-1;                              % followed by [number of LConn, number of RConn]
'Capacitor',[1 2],[1 1],-1;                             % followed by modeling option if any (-1 if none)
%'RLC (Three-Phase)',[1 2 3 4 5 6],[3 3],'ee.passive.rlc_assemblies.rlc.Xabc'; % ABC ports  % modelingOption = get_param(gcbh, 'ComponentPath') 'ee.passive.rlc_assemblies.rlc.abc' get_param(gcb, 'port_option'):'ee.enum.threePhasePort.expanded'
'RLC (Three-Phase)',[1 2 3 4 5 6],[3 3],'ee.enum.threePhasePort.expanded'; % ABC ports  % modelingOption = get_param(gcbh, 'ComponentPath') 'ee.passive.rlc_assemblies.rlc.abc' get_param(gcb, 'port_option'):'ee.enum.threePhasePort.expanded'
'RLC (Three-Phase)',[1 2],[1 1],        'ee.passive.rlc_assemblies.rlc.abc';  % sld ports 
'AC Current Source',[1 2],[1 1],-1;
'AC Voltage Source',[1 2],[1 1],-1;
'DC Voltage Source',[1 2],[1 1],-1;   %10
 ['Voltage' 10 'Source' 10 '(Three-Phase)'],[1 2],[1 1],'ee.sources.voltage.abc';  %SLD
 'Voltage Sensor',[1 3],[1 1],-1;
 'Current Sensor',[1 3],[1 1],-1; 
 [    'Current and Voltage' 10 'Sensor (Three-Phase)'],[1 4],[1 3],'ee.sensors.vi_sensor.abc' ; % SLD %this block has 2 Rconn for outputs meas
 'Phase Splitter',[1 2 3 4],[1,3],-1;  % 4 conns
 ['Circuit Breaker' 10 '(Three-Phase)'],[2 3 4 5 6 7],[3 3],'ee.switches.circuit_breaker.ps.Xabc';
 ['Current' 10 'Source' 10 '(Three-Phase)'],[1 2],[1 1],'ee.sources.current.abc';  %SLD
 ['Coupled Lines' 10 '(Three-Phase)'],[1 2 3 4 5 6],[3 3], 'ee.passive.lines.coupled_lines.Xabc';  % mutual inductance
 ['Nonlinear' 10 'Transformer'],[1 2 3 4],[2 2], -1;   %#19 linear xfo   %21 aout : switch position with mutual inductance
 'Switch',[1 3],[1 1],-1;  %20
 'Voltage Source',[1 2],[1 1],-1;
 'Current Source',[1 2],[1 1],-1;
  ['Circuit Breaker' 10 '(Three-Phase)'],[2 3],[1 1],'ee.switches.circuit_breaker.ps.abc'; %23
  ['Two-Winding' 10 'Transformer' 10 '(Three-Phase)'],[1 2],    [1 1], 'ee.passive.transformers.two_winding_transformer.abc';  %24  NOTE: Y port not available yet
};


SLD_block_index=...
    {[7 11 14 17 23 24 ], ...   %{1} SLD blocks  
    [6 16 18],...        %{2} 3-ph expended blocks  28 aout add 18 mutual inductance!
    [15],...          %{3} phase splitter
    [13 14],...        %{4} current sensors
    [3 4 5 6 7 18 19],...     %{5} rlc branches
    [11],...            % 6: source with RL SLD
    [12],...              %7 voltage sensors
    [-1]};       % 8group of same type, different port option,




nb_blktype=size(blkmask,1);




% find the type of block
qq=[];
for i=1:nb_blk
    flag=0;
    for j=1:nb_blktype
        if isequal(get_param(blk{i,1},'MaskType'),blkmask{j,1})
            if isequal(blkmask{j,4},-1)
                qq=j;
                flag=1;
                break;
            else                
                % cond=isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4});
                % if cond==0
                % try
                %    cond=isequal(get_param(blk{i,1}, 'port_option'),blkmask{j,4});
                % catch
                % end
                % end
                    
                cond1=0; cond2=0;
                try 
                    cond1=isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4}) ;                          
                catch
                   cond1=0;
                end
                try
                    cond2=isequal(get_param(blk{i,1}, 'port_option'),blkmask{j,4}) ; 
                catch
                    cond2=0;
                end
                    
                %if isequal(get_param(blk{i,1}, 'ComponentPath'),blkmask{j,4})
                if cond1==1 | cond2==1
                
                    if isempty(find(SLD_block_index{8}==j))
                            qq=j;
                            flag=1;
                            break;
                    else
                        tmp=get_param(blk{i,1},'Ports');
                        if isequal(tmp(6:7),blkmask{j,3})
                            qq=j;
                            flag=1;
                            break;
                        end
                    end
                end
            end
        end
    end
    if flag==0
        i
        error(['blocktype is not supported yet:' blk{i,1}]);
    end
    blk{i,2}=qq;
    for k=1: length(find(blkmask{qq,2}>0))  % noeud PMIO >0
        blk{i,2+k}=-1;   % initial node number attribution to -1
    end
    if isempty(blk{i,2})
        error('blocktype is null');
    end

end

% reorder to put grounds first
blk1={}; i=1;
blk2={}; j=1;
blk3={}; jj=1;
blk4={}; jjj=1;
cols=size(blk,2);
for k=1:nb_blk
    if (blk{k,2}==1) || (blk{k,2}==2) % ground &  % 3-phase ground
        for m=1:cols
            blk1{i,m}=blk{k,m};
        end
        i=i+1;
    elseif blk{k,2}==19
        for m=1:cols
            blk3{jj,m}=blk{k,m};
        end
        jj=jj+1;
    elseif blk{k,2}==18  %mutual inductance
        for m=1:cols
            blk4{jjj,m}=blk{k,m};
        end
        jjj=jjj+1;

    else
        for m=1:cols
            blk2{j,m}=blk{k,m};
        end
        j=j+1;
    end
end
blk=[blk1 ;blk2; blk3; blk4];
blkgnd=blk1;
%blk=[blk1 ; blk3; blk2];




   

NodeNumber=[];

BlockHandleList=[];
PortHandleList=[];
LineHandleList=[];
NodeNumber=[0];
loop=1;
i=1;  %block indice
jjj=1; %port indice %get ready for non-ground port scan
breakmain=0;
iter=0;
while loop==1
    iter=iter+1;
    if (blk{i,2}==1) || (blk{i,2}==2) % ground &  % 3-phase ground
        PortNum=1;
        NodeNumber=0; 
    else
        % find block and port with -1 node
        % i;   current blk
        % jjj;  current port
if iter==7
    qqqqq=0;
end


        %disp(['i=' num2str(i) ' j=' num2str(jjj)])
        PortNum=blkmask{blk{i,2},2}(jjj);
        NodeNumber=NodeNumber+1;
        typ=get_param(blk{i,1},'BlockType');
        if isequal(typ,'SimscapeBlock')
            blk{i,2+jjj}=NodeNumber;  % if simscape block with -1 : assign the node number now
            % we will find connected blocks after with Portconnectivity
        end
        nexti=i;
        nextjjj=jjj+1;
    end
   
    BlockHandles=[];
    PortHandles=[];
    LineHandles=[];
    SimscapeBlockHandle=[];
    SimscapePortHandle=[];
    %if (blk{i,2}==1) || (blk{i,2}==2) % ground &  % 3-phase ground
        blkh=get(get_param(blk{i,1},'Handle'));
        if (blk{i,2}==1) || (blk{i,2}==2) % ground &  % 3-phase ground
            LineHandles=[LineHandles blkh.LineHandles.LConn ]; % gnd are LConn only (SPECIAL CASE)
        else
            %error('TBDONE')
        end
        %blkh.PortConnectivity; are the cponnections.
        for j=1: length(blkh.PortConnectivity(PortNum).DstBlock)
            % this means that there are some connections
            % Either 'SimscapeBlock','PMIOPort','SubSystem'
            if isequal(get_param(blkh.PortConnectivity(PortNum).DstBlock(j),'Commented'),'off')
                BlockHandles=[BlockHandles blkh.PortConnectivity(PortNum).DstBlock(j)];
                PortHandles=[PortHandles blkh.PortConnectivity(PortNum).DstPort(j)];
                %store the Dest port line handle: part of the current NODE
                tmp=get(blkh.PortConnectivity(PortNum).DstPort(j));
                LineHandles=[LineHandles tmp.Line]; % one line per port
            end
        end

        % BlockHandles=[BlockHandleList];
        % PortHandles=[PortHandleList];
        % LineHandles=[LineHandleList];

        len=length(BlockHandles);
        j=1;
        while j<=len;
            % parcourir BlockHandleList pour trouver d'autres connexions.
            % mais il faut arreter apres 1 passage PMIO->(top-subsystem) ou
            % Subsystem->(inside PMIO)
            % avec
            [new_line_handle, new_block_handle, new_port_handle, flag] = get_more_line4(BlockHandles(j),PortHandles(j),SILENT);
            % if it's a 'SimscapeBlock' , the line has been added previously,
            % then return [][]
            % if it's a 'SubSystem', the outside subsystem line has been added previously,
            % we need to return the PMIO block inside the Subsystem
            if flag==0
            LineHandles=[LineHandles new_line_handle];
            BlockHandles=[BlockHandles new_block_handle];
            PortHandles=[PortHandles new_port_handle];
            disp([blk{i,1} ' Lines--> ' num2str(new_line_handle)]);
            disp([blk{i,1} ' Blocks--> ' num2str(new_block_handle)]);
            len=len+length(new_block_handle);
            else  % if flag==1 -> it is a Simscape block!
                SimscapeBlockHandle=[SimscapeBlockHandle new_block_handle];
                SimscapePortHandle=[SimscapePortHandle new_port_handle];
            end
            j=j+1;
            
        end
    %end
    %PRINT for blk{i}
    if SILENT
    disp(['   ']);
    disp(['block: ' blk{i,1}  ' connections ']);
    end
    for p=1:length(BlockHandles)
        nam=getfullname(BlockHandles(p));
        pget=get(BlockHandles(p));
        allporthan=[pget.PortHandles.LConn pget.PortHandles.RConn];
        ppp=find(PortHandles(p)==allporthan);
        if SILENT
        disp(['blks: ' nam  ' Port: ' num2str(ppp)]);
        end
    end
    for p=1:length(SimscapeBlockHandle)
        nam=getfullname(SimscapeBlockHandle(p));
        pget=get(SimscapeBlockHandle(p));
        allporthan=[pget.PortHandles.LConn pget.PortHandles.RConn];
        ppp=find(SimscapePortHandle(p)==allporthan);
        if SILENT
        disp(['Simscape connected: ' nam  ' Port: ' num2str(ppp)]);
        end
    end
    if SILENT
    disp(['------------------ ']);
    end

    % Add nodes in the blk list
    for p=1:length(SimscapeBlockHandle)
        flagg=0;
        for k=1:nb_blk
            if SimscapeBlockHandle(p) == get_param(blk{k,1},'Handle')
                pget=get(SimscapeBlockHandle(p));
                allporthan=[pget.PortHandles.LConn pget.PortHandles.RConn];
                ppp=find(SimscapePortHandle(p)==allporthan);
                for m=1:length(blkmask{blk{k,2},2})
                   if ppp == blkmask{blk{k,2},2}(m)
                        if blk{k,2+m}<-0.5
                            blk{k,2+m}=NodeNumber;
                            break;
                            flagg=1;
                        end
                    end
                end
                if flagg==1
                    break;
                end
            end
        end
    end
    %increment block-port..
    jjj=jjj+1;
    nbn=length(blkmask{blk{i,2},2});
    if jjj>nbn  % goto next port blk(i)port(jjj)
        jjj=1;
        i=i+1;
        if SILENT
        disp(['block ind ++ i=' num2str(i) ' j=' num2str(jjj)])
        end
    end
    if i>nb_blk
        loop=0;
    else
        while loop==1 && blk{i,2+jjj}>-0.5   % node are -1 for each block , and blk{i,3..4..5}
            nbn=length(blkmask{blk{i,2},2});
            jjj=jjj+1;
            if SILENT
            disp(['node++ i=' num2str(i) ' j=' num2str(jjj)])
            end
            if jjj>nbn  % goto next port blk(i)port(jjj)
                jjj=1;
                i=i+1;
                %disp(['ZZ i=' num2str(i) ' j=' num2str(jjj)])
            end
            if i>nb_blk
                loop=0;
            end
        end
    end




end

%if i>nb_blk
    disp('BEFORE Expension ');
    % add rl of 3 ph source here.
    %TOBE REFINED 
    blk = AddRLsrcNodes(blk,blkmask,SLD_block_index, NodeNumber+4, blkgnd); 
    break_main=1;
    netlistprint(blk);
    disp(' ');
    disp('AFTER SLD Expension ')
    newblk = ExpendSLDNodes(blk,blkmask,SLD_block_index, 1000);  %expand 3-phase blocks
    blk=newblk;
    netlistprint(blk);
    disp(' ');
    disp('AFTER Splitter Expension ')
    blk=MergeSpliterNode(blk,blkmask,SLD_block_index, 1000);
    netlistprint(blk)
    disp(' ');
    disp('AFTER Current sensor node merger ')
    [blk, Vnode,cur_br_all, sensorV,sensorI] = MergeCurrentSensorNodes(blk,SLD_block_index);
    disp('Current sensors and blocks ')
    if SILENT
    for j=1:size(cur_br_all,1)
        disp(['sensor: ' cur_br_all{j,1}]);   
        for k=1:length(cur_br_all{j,2})
            disp(['       right: ' num2str(j) '  ' blk{cur_br_all{j,2}(k),1}]);
        end
        for k=1:length(cur_br_all{j,3})
            disp(['       left: ' num2str(j) '  ' blk{cur_br_all{j,3}(k),1}]);
        end
    end
    end
        netlistprint(blk)
        disp(' ');
        disp('Make SPS netlist ')
        [rlc, source, switches, sensors, rlcnames,sourcenames,source_parameter,switchnames,switchtype,switchVf,switch_parameter] = MakeSPSNetlist(blk,blkmask,1000);
        %netlistprint(blk)
        [yout,y_type,output_sig] = BuildSPSYout(Vnode, switches, cur_br_all,rlcnames,sourcenames,switchnames);
        %y_type=[zeros(1,size(switches,1)) zeros(1,length(Vnode))];  % %%%  added Vnode manually!
        disp('Build MATE interface ')
        [mateblk,source,yout,y_type,sourcenames,nb_nodal_nodes,portType,src_order] = MakeMATEInterface(blk,mateblk_sld,mateblk_1ph,rlc,source,switches,yout,y_type,sourcenames,'findMATEtype');

        disp('Compute ABCD matrices ')
        unit='OMU'; % Units = ohms, mH, and uF
        % %add part of switches on top of source
        % if ~isempty(switches)
        %     swsrc=[switches(:,1:2) ones(size(switches,1),1) zeros(size(switches,1),3)];
        %     source=[swsrc; source];
        % end

        % Open file that contains power_statespace output information
        %fid=fopen('power_circ2ss.net','w');
        %rlc=[rlc(77:end,:); rlc(1:76,:)] ;
         % tmp=1:size(rlc,1);
         % rlc=[rlc tmp'];
         % tmp=ones(size(source,1),1)+1;
         % source=[source tmp];
 
        % rlc=[rlc(end-2:end,:);rlc(1:end-3,:)];
        %[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin]=power_statespace(rlc,switches,source,[],yout,y_type,unit);
        %[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin]=mate_statespace2(rlc,switches,source,[],yout,y_type,unit);
        [A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin,Outputs]=mate_statespace3(rlc,switches,source,[],yout,y_type,unit);
        %Ts=50e-6; should be in the parameter input
        disp('Discretize things... ')
        %nb_nodal_nodes=3;
        %portType=source(end-nb_nodal_nodes+1:end,3);

        %% DISCRETISE METHOD
        if ~isempty(switches)
            [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,switches(:,4)',Ts,nb_nodal_nodes,portType,disc);
        else
            [Adp, B1dp, Cdp, Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p, B2dp]= MATESwitchPermute(A,B,C,D,[],Ts,nb_nodal_nodes,portType,disc);
        end
        statenames=[];   %for later
        u0=zeros(size(source,1),1);
        x0=zeros(size(Adp{1},1),1);  % for later, if required
        disp('Prepare S-function data ')
        MATE=MATE_Sfun_DataFormat(rlc, source, switches, yout,Adp,B1dp,B2dp,Cdp,Ddp, Yp, Cinj_p, Dinj_p, z_p, x_p,mateblk,u0,x0, statenames, rlcnames, sourcenames,switchtype,nb_nodal_nodes,portType, Ts,switchVf);
        MATERT=MATE_Mex_DataFormat(A,B,C,D, switches,switchtype,disc,u0,x0, Ts,MATE.sizes,switchVf);
        disp('Reorder the MATE source from external connexions ')
        try
        RouteMATESources(modele, src_order,'MATE_Equivalent/InputElectricMixer');  %% must change the name if in a subsystem
        catch
            warning('No rewrite order of source done!!!')
        end
        disp('Adding Simulink source in place of Simscape ones ')
        try
        [srcParam,src_Selector]=AddSources(source_parameter,size(switches,1),nb_nodal_nodes);  %srcParam is use directly as a parameter in the model
        catch
            warning('No adds of source done!!!')
        end
        MATE.srcParam=srcParam;
        MATE.src_Selector=src_Selector;
        MATE.source_parameter=source_parameter;
        MATE.switch_parameter=switch_parameter;
        MATE.src_order=src_order;
        MATE.model=modele;
        MATE.Outputs=Outputs;

        disp('Routing Simulink switch signal in place of Simscape ones ') 
        try
        [swParam, masterSwSelOrder]=AddSwitches(modele,switchnames,'InternalSwitchMixer');
        catch
            warning('Error in switch routing... continuing!!!')
        end
        disp('Disabling MATE sources in order to enable simulation from external network')
        MakeMATEInterface([],mateblk_sld,mateblk_1ph,[],[],[],[],[],[],'commentMATE');
        try
            disp('Routing MATE Outputs')   
            SetMATEOutputSignals([model_top '/MATE_Equivalent/MATEOutputSelector'],MATE,sensorV,sensorI,switchnames);
        catch
            warning('Error in output routing... continuing!!!')
        end
    %end


    BlockHandleList=[BlockHandleList BlockHandles];
    PortHandleList=[PortHandleList PortHandles ];
    LineHandleList=[LineHandleList LineHandles];

    % print out

end



