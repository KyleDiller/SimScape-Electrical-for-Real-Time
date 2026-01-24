function [rlc, source, switches, sensors, rlcnames,sourcenames,source_parameter, switchnames, switchtype,switchVf,switch_parameter] = MakeSPSNetlist(blk,maskblk,maxnode)
%[A,B,C,D,states,x0,x0sw,rlsw,u,x,y,freq,Asw,Bsw,Csw,Dsw,Hlin] = power_statespace(rlc,switches,source,line_dist,yout,y_type,unit )
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
rlc=[];
source=[];
sensors=[];
switches=[];
rlcnames={};
sourcenames={};
source_parameter={};
switch_parameter={};
switchnames={};
switchtype=[];
switchVf=[];

no=1;
for i=1:size(blk,1)
    if blk{i,3}>no
        no=blk{i,3};
    end
    if blk{i,4}>no
        no=blk{i,4};
    end
end
MAXNO=maxnode+no;

i=1;
%disp('BUG <= in MAkeSPSNetlist line 24 , Aug 17 2025')
while i<=size(blk,1)
    if blk{i,2}<=2  % forget ground blocks
        i=i+1;
     %%%%%%
     %%% RLC
     %%%%%%%%
    elseif blk{i,2}==3 
        paramNames=[];
        paramNames{1}='R';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        % str='R';
        % try 
        %    val=str2num(get_param(blk{i,1},str));
        % catch
        % end
        % if isempty(val)
        % v=get_param(blk{i,1},'MaskWSVariables');
        % getMaskValue = containers.Map({v.Name}', {v.Value}');
        % val=getMaskValue(str); 
        % end
        rlc=[rlc; blk{i,4} blk{i,3} 0 vals{1}  0     0 ];  % carefull about order node1 node2
        rlcnames{size(rlc,1),1}=blk{i,1};
        rlcnames{size(rlc,1),2}=i;
        i=i+1;
    elseif blk{i,2}==4 
        % str='l';
        % try 
        %    val=str2num(get_param(blk{i,1},str));
        % catch
        % end
        % if isempty(val)
        % v=get_param(blk{i,1},'MaskWSVariables');
        % getMaskValue = containers.Map({v.Name}', {v.Value}');
        % val=getMaskValue(str); 

        paramNames=[];
        paramNames{1}='r';
        paramNames{2}='l';

        vals=ObtainParameterValue(blk{i,1},paramNames);
        rval=vals{1};
        lval=vals{2};
        %disp('Modif MakeSPSNetlist l value 22 sept')
       
        %rlc=[rlc; blk{i,4} blk{i,3} 0  0 val*1e3  0 ];    %L in mH  
        rlc=[rlc; blk{i,4} blk{i,3} 0  rval lval*1e3  0 ];    %L in mH 


        rlcnames{size(rlc,1),1}=blk{i,1};
        rlcnames{size(rlc,1),2}=i;
        i=i+1;
    elseif blk{i,2}==5
        
        % str='c';
        % try 
        %    val=str2num(get_param(blk{i,1},str));
        % catch
        % end
        % if isempty(val)
        % v=get_param(blk{i,1},'MaskWSVariables');
        % getMaskValue = containers.Map({v.Name}', {v.Value}');
        % val=getMaskValue(str); 
        % end
        paramNames=[];
        paramNames{1}='c';
        paramNames{2}='r';
        paramNames{3}='g';

        vals=ObtainParameterValue(blk{i,1},paramNames);
        cval=vals{1};
        rval=vals{2};
        gval=vals{3};

        if gval>0
             error(['MATE blk: ' blk{i,1} 'Capacitor parallel conductance must be equal to 0'])
        end

        %disp('Modif MakeSPSNetlist c value 22 sept')
        %rlc=[rlc; blk{i,4} blk{i,3} 0  0  0   val*1e6 ];    %C in uF 
        rlc=[rlc; blk{i,4} blk{i,3} 0  rval  0   cval*1e6 ];    %C in uF 
        rlcnames{size(rlc,1),1}=blk{i,1};
        rlcnames{size(rlc,1),2}=i;
        i=i+1;

    % elseif blk{i,2}==6
    %     % RLC expended
    %     % val=get_param(blk{i,2},'C');
    %     % rlc=[rlc; blk{i,3} blk{i,4} 0 slResolve( gcb) 0 ];
    %     %TODO
    %     i=i+1;
    elseif blk{i,2}==7 || blk{i,2}==6
        %     RLC SLD      RLC 3ph
        if isequal(get_param(blk{i,1},'MaskType'),maskblk{11,1})
            %this is an added RL from a 3ph source
            paramNames=[];
            paramNames{1}='vline_rms';
            paramNames{2}='shift';
            paramNames{3}='freq';
            paramNames{4}='impedance_option';
            paramNames{5}='R';
            paramNames{6}='L';
            vals=ObtainParameterValue(blk{i,1},paramNames);
            rval=vals{5};
            lval=vals{6};
            cval=0;
            rval=rval*1;   % value in Ohms
            lval=lval*1e3; % value in mH
            cval=0;
            ser=0;
            if isequal(vals{4},0) 
                %no impedance
            else
                if isequal(vals{4},2) 
                    lval=0;
                end
                if isequal(vals{4},3) 
                    rval=0;
                end
                if isequal(vals{4},1) 
                    error(['MATE blk: ' blk{i,1} 'Srouce X/R ratio optnio not supported yet'])
                end


            rlc=[rlc; blk{i,3}   blk{i,4}   ser rval lval cval ];
            rlc=[rlc; blk{i+1,3} blk{i+1,4} ser rval lval cval ];
            rlc=[rlc; blk{i+2,3} blk{i+2,4} ser rval lval cval ];
            rlcnames{size(rlc,1)-2,1}=blk{i,1};
            rlcnames{size(rlc,1)-1,1}=blk{i+1,1};
            rlcnames{size(rlc,1)-0,1}=blk{i+2,1};
            rlcnames{size(rlc,1)-2,2}=i;
            rlcnames{size(rlc,1)-1,2}=i+1;
            rlcnames{size(rlc,1)-0,2}=i+2;
            end
       

        else
            % regular rlc SLD
            rlctype=get_param(blk{i,1}, 'component_structure');
            % rval=slResolve(get_param(blk{i,1},'R'),blk{i,1} ,'expression');
            % lval=slResolve(get_param(blk{i,1},'L'),blk{i,1} ,'expression');
            % cval=slResolve(get_param(blk{i,1},'C'),blk{i,1} ,'expression');
            paramNames=[];
            paramNames{1}='R';
            paramNames{2}='L';
            paramNames{3}='C';
            vals=ObtainParameterValue(blk{i,1},paramNames);
            rval=vals{1};
            lval=vals{2};
            cval=vals{3};


            if isequal(rlctype,'ee.enum.rlc.structure.R')
                ser=0; lval=0; cval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.L')
                ser=0; rval=0; cval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.C')
                ser=0; lval=0; rval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRL')
                ser=0; cval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRC')
                ser=0; lval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesLC')
                ser=0; rval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRLC')
                ser=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRL')
                ser=1; cval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRC')
                ser=1; lval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelLC')
                ser=1; rval=0;
            elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRLC')
                ser=1;
            else
            end
            rval=rval*1;   % value in Ohms
            lval=lval*1e3; % value in mH
            cval=cval*1e6; % value in uF
            rlc=[rlc; blk{i,3}   blk{i,4}   ser rval lval cval ];
            rlc=[rlc; blk{i+1,3} blk{i+1,4} ser rval lval cval ];
            rlc=[rlc; blk{i+2,3} blk{i+2,4} ser rval lval cval ];
            rlcnames{size(rlc,1)-2,1}=blk{i,1};
            rlcnames{size(rlc,1)-1,1}=blk{i+1,1};
            rlcnames{size(rlc,1)-0,1}=blk{i+2,1};
            rlcnames{size(rlc,1)-2,2}=i;
            rlcnames{size(rlc,1)-1,2}=i+1;
            rlcnames{size(rlc,1)-0,2}=i+2;
        end
        i=i+3;
        %%%%%%%%%%%%%%%
        %%% TRANSFO et MUTUAL INDUCTANCE
        %%%%%%%%%%%%%%%
    elseif blk{i,2}==19  % xfo
        
        paramNames=[];
        paramNames{1}='Nw';
        paramNames{2}='Nw2';
        paramNames{3}='R_1';
        paramNames{4}='L_1';
        paramNames{5}='R_2';
        paramNames{6}='L_2';
        paramNames{7}='R_m';
        paramNames{8}='parameterization_option';
        paramNames{9}='L';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        if ~isequal(vals{8},1)
            error(['MATE blk: ' blk{i,1} 'Only linear xfo is suported now'])
        end
        
        rlc=[rlc; blk{i,3}             blk{i,4}  2 vals{3} vals{4}*1e3  vals{1}];
        rlc=[rlc; blk{i,5}             blk{i,6}  2 vals{5} vals{6}*1e3  vals{2}];
        %rlc=[rlc; blk{i,5}+MAXNO      blk{i,4}  1 vals{7} vals{9}*1e3     0   ];  % internal node
        rlc=[rlc; MAXNO      blk{i,4}  1 vals{7} vals{9}*1e3     0   ];
        MAXNO=MAXNO+1;
        rlcnames{size(rlc,1)-2,1}=[blk{i,1} '_PRIM'];
        rlcnames{size(rlc,1)-1,1}=[blk{i,1} '_SEC'];
        rlcnames{size(rlc,1)-0,1}=[blk{i,1} '_MAG'];
        rlcnames{size(rlc,1)-2,2}=i;
        rlcnames{size(rlc,1)-1,2}=i;
        rlcnames{size(rlc,1)-0,2}=i;
        i=i+1;
        
     elseif blk{i,2}==18 % mutual inductance
        paramNames=[];
        paramNames{1}='prm';
        paramNames{2}='La';
        paramNames{3}='Lb';
        paramNames{4}='Lc';
        paramNames{5}='Lmab';
        paramNames{6}='Lmac';
        paramNames{7}='Lmbc';
        paramNames{8}='Ra';
        paramNames{9}='Rb';
        paramNames{10}='Rc';
        paramNames{11}='Rm';
        paramNames{12}='Ga';
        paramNames{13}='Gb';
        paramNames{14}='Gc';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        if ~isequal(vals{1},2)
            error(['MATE blk: ' blk{i,1} 'Only Generalized inductance supported is suported now'])
        end
   %  6.0000         1    4.0000    1.0000    1.0000    0    1.0000 aa
   %  5.0000         2    4.0000    1.0000    1.0000    0    2.0000 bb
   %  4.0000         3    4.0000    1.0000    1.0000    0    3.0000 cc
   % 12.0000 uncon   1         0    0.9000    0.9000    0    4.0000 ab
   % 13.0000 uncon   1         0    0.9000    0.8000    0    5.0000 ac
   % 14.0000 uncon   2         0    0.9000    0.7000    0    6.0000 bc  
   % 1 1 2 is not an error.

        rlc=[rlc; blk{i,3}             blk{i,4}    4 vals{8}  vals{2}*1e3  0];
        rlc=[rlc; blk{i+1,3}           blk{i+1,4}  4 vals{9}  vals{3}*1e3  0];
        rlc=[rlc; blk{i+2,3}           blk{i+2,4}  4 vals{10} vals{4}*1e3  0];
        % rlc=[rlc; blk{i,3}+MAXNO       blk{i,4}    0 vals{11} vals{5}*1e3  0];  % internal node
        % rlc=[rlc; blk{i,3}+MAXNO+1     blk{i,4}    0 vals{11} vals{6}*1e3  0];  % internal node
        % rlc=[rlc; blk{i,3}+MAXNO+2     blk{i+1,4}  0 vals{11} vals{7}*1e3  0];  % internal node
        rlc=[rlc; MAXNO       blk{i,4}    0 vals{11} vals{5}*1e3  0];  % internal node
        rlc=[rlc; MAXNO+1     blk{i,4}    0 vals{11} vals{6}*1e3  0];  % internal node
        rlc=[rlc; MAXNO+2     blk{i+1,4}  0 vals{11} vals{7}*1e3  0];  % internal node
        MAXNO=MAXNO+3;
        w=size(rlc,1);
        rlcnames{w-5,1}=[blk{i,1} '_Maa'];
        rlcnames{w-4,1}=[blk{i,1} '_Mbb'];
        rlcnames{w-3,1}=[blk{i,1} '_Mcc'];
        rlcnames{w-2,1}=[blk{i,1} '_Mab'];
        rlcnames{w-1,1}=[blk{i,1} '_Mac'];
        rlcnames{w-0,1}=[blk{i,1} '_Mbc'];
        rlcnames{w-5,2}=i;
        rlcnames{w-4,2}=i+1;
        rlcnames{w-3,2}=i+2;
        rlcnames{w-2,2}=-i;
        rlcnames{w-1,2}=-i;
        rlcnames{w-0,2}=-i;

        disp(['Mutual inductance ' blk{i,1} ' : neglecting shunt conductance']);

        i=i+3;
 % 28 aout: block is triple in blk{} so we treat the 1st occurance only

  %Two-Winding' 10 'Transformer (Three-Phase) SLD 'ee.passive.transformers.two_winding_transformer.abc';  %24
  elseif blk{i,2}==24
        paramNames=[];
        paramNames{1}='SRated';
        paramNames{2}='FRated';
        paramNames{3}='Winding1Connection'; 
                                            % ee.enum.windingconnection.Y
                                            % ee.enum.windingconnection.Yn
                                            % ee.enum.windingconnection.Yg
                                            % ee.enum.windingconnection.delta1
                                            % ee.enum.windingconnection.delta11
        paramNames{4}='VRated1';
        paramNames{5}='Winding2Connection';
        paramNames{6}='VRated2';
        paramNames{7}='CoreType'; %ee.enum.coretype.threelimb ee.enum.coretype.fivelimb
        paramNames{8}='pu_Rw1';
        paramNames{9}='pu_Rw2';
        paramNames{10}='pu_Xl1';
        paramNames{11}='pu_Xl2';
        paramNames{12}='magnetizing_resistance_option'; %ee.enum.transformer_magnetizingResistance.exclude ee.enum.transformer_magnetizingResistance.include
        paramNames{13}='pu_Rm';
        paramNames{14}='magnetizing_reactance_option'; %ee.enum.transformer_magnetizingReactance.exclude ee.enum.transformer_magnetizingReactance.include
        paramNames{15}='pu_Xm';
        paramNames{16}='saturation_option'; %ee.enum.transformer_saturation.exclude ee.enum.transformer_saturation.include

        vals=ObtainParameterValue(blk{i,1},paramNames);
        if isequal(vals{3},'Yn')  | isequal(vals{5},'Yn')
            error(['MATE blk: ' blk{i,1} 'Xfo neutral connecton not available. Please use single-phase xfo to model this.'])
        end
        tmp1=char(vals{3}); tmp2=char(vals{5});
        if isequal(tmp1(1),'d')  | isequal(tmp2(1),'d')  % for 'delta'
            warning(['MATE blk: ' blk{i,1} 'Xfo Delta connection not available. Please use single-phase xfo to model this.'])
        end
        if ~isequal(vals{7},'fivelimb')
            error(['MATE blk: ' blk{i,1} '5-limb 3-phase xfo supported only'])
        end
        if ~isequal(vals{12},'include')
            error(['MATE blk: ' blk{i,1} '3-phase xfo: must include magnetizing resistance'])
        end
        if ~isequal(vals{14},'include')
            error(['MATE blk: ' blk{i,1} '3-phase xfo: must include magnetizing inductance'])
        end
        if ~isequal(vals{16},'exclude')
            error(['MATE blk: ' blk{i,1} '3-phase xfo: saturation option not supported (yet)'])
        end
        if isequal(vals{3},'Y') 
            noY1=MAXNO;MAXNO=MAXNO+1;
        elseif isequal(vals{3},'Yg') 
            noY1=0;
        else
        end
        if isequal(vals{5},'Y') 
            noY2=MAXNO;MAXNO=MAXNO+1;
        elseif isequal(vals{5},'Yg') 
            noY2=0;
        else
        end
        %pu->SI conversion
        Sbase=vals{1}/3;
        wbase=2*pi*vals{2};
        if isequal(vals{3},'Y')  | isequal(vals{3},'Yg') || isequal(vals{3},'Yn')
            Vbase=vals{4}/sqrt(3);
        else
            Vbase=vals{4};
        end
        Ibase=Sbase/Vbase;
        Zbase=Vbase/Ibase;
        Lbase=Zbase/wbase;

        if isequal(vals{5},'Y')  | isequal(vals{5},'Yg') || isequal(vals{5},'Yn')
            Vbase2=vals{6}/sqrt(3);
        else
            Vbase2=vals{6};
        end
        Ibase2=Sbase/Vbase2;
        Zbase2=Vbase2/Ibase2;
        Lbase2=Zbase2/wbase;

        if isequal(vals{3},'Y')  | isequal(vals{3},'Yg') || isequal(vals{3},'Yn')
                connP=1;
            elseif isequal(vals{3},'delta1')
                connP=2;
            elseif isequal(vals{3},'delta11')
                connP=3;
        else
            error('xfo3:1')
        end
        if isequal(vals{5},'Y')  | isequal(vals{5},'Yg') || isequal(vals{5},'Yn')
                connS=1;
            elseif isequal(vals{5},'delta1')
                connS=2;
            elseif isequal(vals{5},'delta11')
                connS=3;
        else
            error('xfo3:2')
        end

        if connP==1
            rlc=[rlc; blk{i,3}     noY1        2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        elseif connP==2
            rlc=[rlc; blk{i,3}     blk{i+1,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        else
            rlc=[rlc; blk{i,3}     blk{i+2,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        end
        if connS==1
            rlc=[rlc; blk{i,4}     noY2          2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        elseif connS==2
            rlc=[rlc; blk{i,4}     blk{i+1,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        else
            rlc=[rlc; blk{i,4}     blk{i+2,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        end
        rlc=[rlc;  MAXNO       rlc(end-1,2)  1 Zbase*vals{13}  Lbase*vals{15}*1e3   0];  % internal node
        MAXNO=MAXNO+1;

        rlcnames{size(rlc,1)-2,1}=[blk{i,1} '_PRIM1'];
        rlcnames{size(rlc,1)-1,1}=[blk{i,1} '_SEC1'];
        rlcnames{size(rlc,1)-0,1}=[blk{i,1} '_MAG1'];
        rlcnames{size(rlc,1)-2,2}=i;
        rlcnames{size(rlc,1)-1,2}=i;
        rlcnames{size(rlc,1)-0,2}=i;


        % rlc=[rlc; blk{i+1,3}     noY1  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        % rlc=[rlc; blk{i+1,4}     noY2  2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        % rlc=[rlc;  MAXNO         noY1  1 Zbase*vals{13}  Lbase*vals{15}*1e3   0];  % internal node
        % MAXNO=MAXNO+1;
        if connP==1
            rlc=[rlc; blk{i+1,3}     noY1        2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        elseif connP==2
            rlc=[rlc; blk{i+1,3}     blk{i+2,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        else
            rlc=[rlc; blk{i+1,3}     blk{i,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        end
        if connS==1
            rlc=[rlc; blk{i+1,4}     noY2          2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        elseif connS==2
            rlc=[rlc; blk{i+1,4}     blk{i+2,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        else
            rlc=[rlc; blk{i+1,4}     blk{i,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        end
        rlc=[rlc;  MAXNO       rlc(end-1,2)      1 Zbase*vals{13}  Lbase*vals{15}*1e3   0];  % internal node
        MAXNO=MAXNO+1;

        rlcnames{size(rlc,1)-2,1}=[blk{i,1+1} '_PRIM2'];
        rlcnames{size(rlc,1)-1,1}=[blk{i,1+1} '_SEC2'];
        rlcnames{size(rlc,1)-0,1}=[blk{i,1+1} '_MAG2'];
        rlcnames{size(rlc,1)-2,2}=i;
        rlcnames{size(rlc,1)-1,2}=i;
        rlcnames{size(rlc,1)-0,2}=i;
        % 
        % rlc=[rlc; blk{i+2,3}     noY1  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        % rlc=[rlc; blk{i+2,4}     noY2  2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        % rlc=[rlc;  MAXNO         noY1  1 Zbase*vals{13}  Lbase*vals{15}*1e3   0];  % internal node
        % MAXNO=MAXNO+1;

        if connP==1
            rlc=[rlc; blk{i+2,3}     noY1        2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        elseif connP==2
            rlc=[rlc; blk{i+2,3}     blk{i,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        else
            rlc=[rlc; blk{i+2,3}     blk{i+1,3}  2 Zbase*vals{8}   Lbase*vals{10}*1e3   Vbase];
        end
        if connS==1
            rlc=[rlc; blk{i+2,4}     noY2          2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        elseif connS==2
            rlc=[rlc; blk{i+2,4}     blk{i,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        else
            rlc=[rlc; blk{i+2,4}     blk{i+1,4}    2 Zbase2*vals{9}  Lbase2*vals{11}*1e3  Vbase2];
        end
        rlc=[rlc;  MAXNO       rlc(end-1,2)      1 Zbase*vals{13}  Lbase*vals{15}*1e3   0];  % internal node
        MAXNO=MAXNO+1;

        rlcnames{size(rlc,1)-2,1}=[blk{i,1+2} '_PRIM3'];
        rlcnames{size(rlc,1)-1,1}=[blk{i,1+2} '_SEC3'];
        rlcnames{size(rlc,1)-0,1}=[blk{i,1+2} '_MAG3'];
        rlcnames{size(rlc,1)-2,2}=i;
        rlcnames{size(rlc,1)-1,2}=i;
        rlcnames{size(rlc,1)-0,2}=i;

        i=i+3;

        %%%%%%%%%%%%
        %%% SOURCE
        %%%%%%%%%%%%
    elseif blk{i,2}==8 %AC Current Source
        paramNames=[];
        paramNames{1}='amp';
        paramNames{2}='shift';
        paramNames{3}='frequency';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        source=[source ; blk{i,3} blk{i,4} 1 vals{1} vals{2} vals{3}];
        sourcenames{size(source,1),1}=blk{i,1};
        sourcenames{size(source,1),2}=i;
        len=size(source,1);
        source_parameter{len,1}=blk{i,2};
        source_parameter{len,2}=[vals{1} vals{2} vals{3}];
        i=i+1;
    elseif blk{i,2}==9 %AC voltage Source
        paramNames=[];
        paramNames{1}='amp';
        paramNames{2}='shift';
        paramNames{3}='frequency';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        source=[source ; blk{i,3} blk{i,4} 0 vals{1} vals{2} vals{3}];
        sourcenames{size(source,1),1}=blk{i,1};
        sourcenames{size(source,1),2}=i;
        len=size(source,1);
        source_parameter{len,1}=blk{i,2};
        source_parameter{len,2}=[vals{1} vals{2} vals{3}];
        i=i+1;
    elseif blk{i,2}==10 %DC voltage Source
        paramNames=[];
        paramNames{1}='amp';
        paramNames{2}='shift';
        paramNames{3}='frequency';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        source=[source ; blk{i,3} blk{i,4} 0 vals{1} vals{2} vals{3}];
        sourcenames{size(source,1),1}=blk{i,1};
        sourcenames{size(source,1),2}=i;
        len=size(source,1);
        source_parameter{len,1}=blk{i,2};
        source_parameter{len,2}=[vals{1} vals{2} vals{3}];
        i=i+1;
    elseif blk{i,2}==11 % 3 phase voltage Source SLD
        paramNames=[];
        paramNames{1}='vline_rms';
        paramNames{2}='shift';
        paramNames{3}='freq';
        paramNames{4}='impedance_option';
        paramNames{5}='R';
        paramNames{6}='L';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        if vals{4}==0
            rr=0;
            ll=0;
        elseif vals{4}==2
            rr=vals{5};  % not sure why this is a number and not a string...
            ll=0;
        elseif vals{4}==4
            rr=vals{5};
            ll=vals{6};
        else
            error(['MakeNetlist: the selected source impedance option is not supported. Block: ' blk{i,1}]);
        end
        if isempty(source)
            allnodes=unique([rlc(:,1)' rlc(:,2)']);
        else
            allnodes=unique([rlc(:,1)' rlc(:,2)' source(:,1)' source(:,2)']);
        end
        newnode=max(allnodes)+1;
        % if ll==0 && rr==0  % no more condition because RL is any is
        % treated as a separated RL added blk{}
        x=-90;
            source=[source ; blk{i,4}   blk{i,3}   0 vals{1}*sqrt(2/3) vals{2}     vals{3}];
            source=[source ; blk{i+1,4} blk{i+1,3} 0 vals{1}*sqrt(2/3) vals{2}-120 vals{3}];
            source=[source ; blk{i+2,4} blk{i+2,3} 0 vals{1}*sqrt(2/3) vals{2}+120 vals{3}];  %%%check out the node 1 node 2 order!!
            sourcenames{size(source,1)-2,1}=blk{i,1};
            sourcenames{size(source,1)-1,1}=blk{i+1,1};
            sourcenames{size(source,1)-0,1}=blk{i+2,1};
            sourcenames{size(source,1)-2,2}=i;
            sourcenames{size(source,1)-1,2}=i+1;
            sourcenames{size(source,1)-0,2}=i+2;
        len=size(source,1);
        source_parameter{len-2,1}=blk{i,2};
        source_parameter{len-2,2}=[vals{1}*sqrt(2/3) vals{2}     vals{3}];
        source_parameter{len-1,1}=blk{i,2};
        source_parameter{len-1,2}=[vals{1}*sqrt(2/3) vals{2}-120 vals{3}];
        source_parameter{len-0,1}=blk{i,2};
        source_parameter{len-0,2}=[vals{1}*sqrt(2/3) vals{2}+120 vals{3}];

        % else
        %     % add an intermediate node to insert the RL
        %     source=[source ; blk{i,3}   newnode   0 vals{1}*sqrt(2/3) vals{2}     vals{3}];
        %     source=[source ; blk{i+1,3} newnode+1 0 vals{1}*sqrt(2/3) vals{2}+120 vals{3}];
        %     source=[source ; blk{i+2,3} newnode+2 0 vals{1}*sqrt(2/3) vals{2}-120 vals{3}];
        %     sourcenames{size(source,1)-2,1}=blk{i,1};
        %     sourcenames{size(source,1)-1,1}=blk{i+1,1};
        %     sourcenames{size(source,1)-0,1}=blk{i+2,1};
        %     sourcenames{size(source,1)-2,2}=i;
        %     sourcenames{size(source,1)-1,2}=i+1;
        %     sourcenames{size(source,1)-0,2}=i+2;
        %     rlc=[rlc; newnode   blk{i,4}   0 rr ll 0 ];
        %     rlc=[rlc; newnode+1 blk{i+1,4} 0 rr ll 0 ];
        %     rlc=[rlc; newnode+2 blk{i+2,4} 0 rr ll 0 ]; % vadd source impedance
        % rlcnames{size(rlc,1)-2,1}=blk{i,1};
        % rlcnames{size(rlc,1)-1,1}=blk{i+1,1};
        % rlcnames{size(rlc,1)-0,1}=blk{i+2,1};
        % rlcnames{size(rlc,1)-2,2}=i;
        % rlcnames{size(rlc,1)-1,2}=i+1;
        % rlcnames{size(rlc,1)-0,2}=i+2;
        %end
        i=i+3;
      elseif blk{i,2}==17 % 3 phase current Source SLD
        paramNames=[];
        paramNames{1}='iphase_rms';
        paramNames{2}='shift';
        paramNames{3}='freq';

        vals=ObtainParameterValue(blk{i,1},paramNames);

            source=[source ; blk{i,3}   blk{i,4}   1 vals{1}*sqrt(2/3) vals{2}     vals{3}];
            source=[source ; blk{i+1,3} blk{i+1,4} 1 vals{1}*sqrt(2/3) vals{2}-120 vals{3}];
            source=[source ; blk{i+2,3} blk{i+2,4} 1 vals{1}*sqrt(2/3) vals{2}+120 vals{3}];
            sourcenames{size(source,1)-2,1}=blk{i,1};
            sourcenames{size(source,1)-1,1}=blk{i+1,1};
            sourcenames{size(source,1)-0,1}=blk{i+2,1};
            sourcenames{size(source,1)-2,2}=i;
            sourcenames{size(source,1)-1,2}=i+1;
            sourcenames{size(source,1)-0,2}=i+2;
        len=size(source,1);
        source_parameter{len-2,1}=blk{i,2};
        source_parameter{len-2,2}=[vals{1}*sqrt(2/3) vals{2}     vals{3}];
        source_parameter{len-1,1}=blk{i,2};
        source_parameter{len-1,2}=[vals{1}*sqrt(2/3) vals{2}-120 vals{3}];
        source_parameter{len-0,1}=blk{i,2};
        source_parameter{len-0,2}=[vals{1}*sqrt(2/3) vals{2}+120 vals{3}];
        i=i+3;
  %%%%%%%%%%%%%%
  %%%% switches
  %%%%%%%%%%%%%%%%
  %[ node1, node2, status, R, L/Xl, no_I , no_U ] 
    elseif blk{i,2}==16 || blk{i,2}==23 % 3phase breaker
        paramNames=[];
        paramNames{1}='R_closed';
        paramNames{2}='G_open';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        nsw=size(switches,1);
        switches=[switches ; blk{i,3}   blk{i,4}   0 vals{1} 0 nsw+1 nsw+1];
        switches=[switches ; blk{i+1,3} blk{i+1,4} 0 vals{1} 0 nsw+2 nsw+2];
        switches=[switches ; blk{i+2,3} blk{i+2,4} 0 vals{1} 0 nsw+3 nsw+3];
        
        rlc=[rlc; blk{i,3}   blk{i,4}   0 1/vals{2} 0 0 ];
        rlc=[rlc; blk{i+1,3} blk{i+1,4} 0 1/vals{2} 0 0 ];
        rlc=[rlc; blk{i+2,3} blk{i+2,4} 0 1/vals{2} 0 0 ]; % very important to add this Ropen in parallel to switch ottherwise it does not compute...
        rlcnames{size(rlc,1)-2,1}=blk{i,1};
        rlcnames{size(rlc,1)-1,1}=blk{i+1,1};
        rlcnames{size(rlc,1)-0,1}=blk{i+2,1};
        rlcnames{size(rlc,1)-2,2}=i;
        rlcnames{size(rlc,1)-1,2}=i+1;
        rlcnames{size(rlc,1)-0,2}=i+2;
        switchnames{size(switches,1)-2,1}=blk{i,1};
        switchnames{size(switches,1)-1,1}=blk{i,1};
        switchnames{size(switches,1)-0,1}=blk{i,1};  % put all reference to the first block
        switchnames{size(switches,1)-2,2}=i;
        switchnames{size(switches,1)-1,2}=i+1;
        switchnames{size(switches,1)-0,2}=i+2;
        switchtype=[switchtype 2 2 2];
        switchVf=[switchVf 0 0 0];
        len=size(switches,1);
        switch_parameter{len-2,1}=blk{i,2};
        switch_parameter{len-2,2}=[vals{1} 0];  % Ron Vf
        switch_parameter{len-1,1}=blk{i,2};
        switch_parameter{len-1,2}=[vals{1} 0];
        switch_parameter{len-0,1}=blk{i,2};
        switch_parameter{len-0,2}=[vals{1} 0];
   
        i=i+3;
    elseif blk{i,2}==20 % 1 phase switch
        paramNames=[];
        paramNames{1}='R_closed';
        paramNames{2}='G_open';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        nsw=size(switches,1);
        switches=[switches ; blk{i,3}   blk{i,4}   0 vals{1} 0 nsw+1 nsw+1];
        rlc=[rlc; blk{i,3}   blk{i,4}   0 1/vals{2} 0 0 ]; % very important to add this Ropen in parallel to switch ottherwise it does not compute..
        rlcnames{size(rlc,1),1}=blk{i,1};
        rlcnames{size(rlc,1),2}=i;

        switchnames{size(switches,1),1}=blk{i,1};
        switchnames{size(switches,1),2}=i;
        switchtype=[switchtype 1];
        switchVf=[switchVf 0];
        len=size(switches,1);
        switch_parameter{len,1}=blk{i,2};
        switch_parameter{len,2}=[vals{1} 0];  % Ron Vf
        i=i+1;
    else
        i=i+1;
    end
end
% check for short-ciruit of rlc
for i=1:size(rlc,1)
    if (rlc(i,1)==rlc(i,2))
        error(['MAKESPSNETLIST 168: short circuit in block: ' blk{i,1}])
    end
end

%add part of switches on top of source
if ~isempty(switches)
    swsrc=[switches(:,1:2) ones(size(switches,1),1) zeros(size(switches,1),3)];
    source=[swsrc; source];
end
% add sourcenames for these
% nb_switch=size(switches,1);
% for i=size(sourcenames,1):-1:1
%     sourcenames{i+nb_switch,1}=sourcenames{i,1};
%     sourcenames{i+nb_switch,2}=sourcenames{i,2};
% end
% for i=1:nb_switch
%     sourcenames{i,1}=['sw' num2str(i)];
%     sourcenames{i,2}=[-i];
% end

sourcenames=[switchnames ; sourcenames];

%record the real internal source params  

% for i=nb_switch+1:size(sourcenames,1)  % MATE source not included yet
% 
% end




end

    
function vals=ObtainParameterValue(blk,paramNames)

vals=cell(length(paramNames),1);

for i=1:length(paramNames)
        str=paramNames{i};
        try 
        v=get_param(blk,'MaskWSVariables');
        getMaskValue = containers.Map({v.Name}', {v.Value}');
        val=getMaskValue(str); 
        vals{i}=val;
        if isempty(val)
            xxx=get_param(blk,'MaskNames');
            yyy=get_param(blk,'MaskValues');
            zzz=containers.Map(xxx, yyy);
            val=zzz(str);
            vals{i}=evalin('base',val);
        end
        catch
            error(['3: not able to obtain params of block' blk])
        end
end
end
