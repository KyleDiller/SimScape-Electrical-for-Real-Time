function [mateblk,source,yout,y_type, sourcenames,nb_nodal_nodes,portType,src_order] = MakeMATEInterface(blk,mateblk_sld,mateblk_1ph,rlc,source,switches,yout,y_type,sourcenames, step)
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

% wow! since we add the nodal source and output at the end, they appear at
% the end fo in-out vectors!
%nb_nodal_nodes=3;
%portType=source(end-nb_nodal_nodes+1:end,3)
if isempty(mateblk_1ph)
    mateblk=[mateblk_sld];
elseif isempty(mateblk_sld)
    mateblk=[mateblk_1ph];
else
    mateblk=[mateblk_sld; mateblk_1ph];
end
ind_sld=size(mateblk_sld,1);

youtV=[];
y_typeV=[];
portTypeV=[];
youtI=[];
y_typeI=[];
portTypeI=[];
sourceV=[];
sourcenamesV=[];
sourceI=[];
sourcenamesI=[];

if isequal(step,'findMATEtype')
    nb_nodal_nodes=0;
    portType=[];
    for m=1:size(mateblk,1)
        mateblk{m,2}=get_param(mateblk{m,1},'bn');
        flag=0;
        for j=1:size(blk,1)
            if isequal([mateblk{m,1} '/MATEshunt'],blk{j,1})
                if (blk{j,3}==0)
                    node=blk{j,4};
                    index=[j 4];
                else
                    node=blk{j,3}; % node is the one of the MATE interface
                    index=[j 3];
                end

                flag=1;
            end
        end
        if flag<1
            error('MakeMATEInterface 1')
        end
        flag=0;
        for j=1:size(blk,1)
            if isequal([mateblk{m,1} '/MATEserie'],blk{j,1})
                if (blk{j,3}==node)
                    nodenot=blk{j,4};
                else
                    nodenot=blk{j,3}; % nodenot is the node that must NOT be considered in MATe type search
                end
                flag=flag+1;
            end
        end
        if flag<1
            error('MakeMATEInterface 2')
        end

        % trouver chemin de node vers le ground... si oui, I-type
        thres=500; % valeur critique Thevenin/Norton pour le cas Ronly , en ohms
        type=1; % voltage or capacitor path -> path to ground possible
        type=0; % current or inductor path -> not path to ground
        allnodes=[];
        for i=1:size(rlc,1)
            if rlc(i,5)==0 && rlc(i,6)==0
                if rlc(i,4)>thres
                    type=0;
                else
                    type=1; % cas Ronly
                end
            elseif rlc(i,3)==1  % parallel branch
                if rlc(i,6)==0
                    type=0;  % cas parallel sans C
                else
                    type=1;
                end
            elseif rlc(i,3)==0  %series branches
                if rlc(i,5)==0
                    type=1;  % cas serie sans L
                else
                    type=0;
                end
            else
                type=0;
            end
            allnodes=[allnodes; rlc(i,1:2) type -1 ];
        end
        for i=1:size(source,1)
            if source(i,3)==0
                type=1;  % means a voltage source do not block a path to ground
            else
                type=0;
            end
            allnodes=[allnodes; source(i,1:2) type -1 ];
        end

        % mask as read the rlc source with node 'nodenot'
        % because we want to serach ground path in the direction of the
        % MATE group
        for j=1:size(allnodes,1)
            if allnodes(j,1)==nodenot || allnodes(j,2)==nodenot
                allnodes(j,4)=1;
            end
        end


        i=1;
        nodeX=node;
        while i<=length(nodeX)
            for j=1:size(allnodes,1)
                if allnodes(j,1)==nodeX(i) || allnodes(j,2)==nodeX(i)
                    if allnodes(j,3)==1 && allnodes(j,4)<0
                        if allnodes(j,1)==nodeX(i)
                            nodeX=[nodeX allnodes(j,2)]; % add the other branch node
                        else
                            nodeX=[nodeX allnodes(j,1)]; % add the  branch node
                        end
                        nodeX=unique(nodeX);
                        allnodes(j,4)=1; % mask as read
                    end
                end
            end
            i=i+1;
        end
        if ~isempty(find(nodeX==0))
            % path to ground found
            MATEtype='I';
        else
            MATEtype='V';
        end
        mateblk{m,3}=MATEtype;
        if MATEtype=='V'
            set_param([mateblk{m,1} '/MATE_V'],'Commented','off');
            set_param([mateblk{m,1} '/MATE_I'],'Commented','on');
        else
            set_param([mateblk{m,1} '/MATE_V'],'Commented','on');
            set_param([mateblk{m,1} '/MATE_I'],'Commented','off');
        end

        % ajouter les sorties MATE
        % il faut trouver les branches Ibn dans rlc et In dans source
        rlcbr=zeros(size(rlc,1),2);
        srcbr=zeros(size(source,1),2);
        % find branches that are connected to node and not connected to
        % nodenot (nodenot is connected to the MATEserie branch)
        if MATEtype=='V'  % than we must find the output current
            for i=1:size(rlc,1)
                if rlc(i,1)==node || rlc(i,2)==node
                    if rlc(i,1)~=nodenot && rlc(i,2)~=nodenot
                        rlcbr(i,1)=1;
                        if rlc(i,2)==node
                            rlcbr(i,2)='-';
                        else
                            rlcbr(i,2)='+';
                        end
                    end
                end
            end
            for i=1:size(source,1)
                if source(i,1)==node || source(i,2)==node
                    if source(i,1)~=nodenot && source(i,2)~=nodenot
                        srcbr(i,1)=1;
                        if source(i,2)==node
                            srcbr(i,2)='-';
                        else
                            srcbr(i,2)='+';
                        end
                    end
                end
            end

            % built the yout expression
            str=[];strb=[]; strc=[];
            pos=find(rlcbr==1);
            for i=1:length(pos)

                    str =[str  rlcbr(pos(i),2)];
                    if m<=ind_sld
                        strb=[strb rlcbr(pos(i),2)];
                        strc=[strc rlcbr(pos(i),2)];
                    end

                str= [str  'I_b' num2str(pos(i)+0)];
                if m<=ind_sld
                    strb=[strb 'I_b' num2str(pos(i)+1)];
                    strc=[strc 'I_b' num2str(pos(i)+2)];
                end
            end

            pos2=find(srcbr==1);
            if ~isempty(str) && ~isempty(pos2)
                str =[str  srcbr(pos2(i),2)];
                if m<=ind_sld
                    strb=[strb srcbr(pos2(i),2)];
                    strc=[strc srcbr(pos2(i),2)];
                end
            end
            for i=1:length(pos2)
                if i>1
                    str =[str  srcbr(pos2(i),2)];
                    if m<=ind_sld
                        strb=[strb srcbr(pos2(i),2)];
                        strc=[strc srcbr(pos2(i),2)];
                    end
                end
                str =[str  'I' num2str(pos2(i)+0)];
                if m<=ind_sld
                    strb=[strb 'I' num2str(pos2(i)+1)];
                    strc=[strc 'I' num2str(pos2(i)+2)];
                end
            end
            %%% build for 3ph
            if m<=ind_sld
                % if isempty(yout)
                %     yout=char(str, strb, strc);
                % else
                %     yout=char(yout,str, strb, strc);
                % end
                % y_type=[y_type 1 1 1];
                % nb_nodal_nodes=nb_nodal_nodes+3;
                % portType=[portType 0 0 0];

                if isempty(youtV)
                    youtV=char(str, strb, strc);
                else
                    youtV=char(youtV,str, strb, strc);
                end
                y_typeV=[y_typeV 1 1 1];
                portTypeV=[portTypeV 0 0 0];
                nb_nodal_nodes=nb_nodal_nodes+3;
            else
                % if isempty(yout)
                %     yout=char(str);
                % else
                %     yout=char(yout,str);
                % end
                % y_type=[y_type 1];
                % nb_nodal_nodes=nb_nodal_nodes+1;
                % portType=[portType 0];
                if isempty(youtV)
                    youtV=char(str);
                else
                    youtV=char(youtV,str);
                end
                y_typeV=[y_typeV 1];
                portTypeV=[portTypeV 0];
                nb_nodal_nodes=nb_nodal_nodes+1;
            end

        end

        if MATEtype=='I'  % then we must find the node voltage: EASY
            str=['U_n' num2str(node) '_0'];
            % finf voltage of B and C node
            if m<=ind_sld
                nodeb=blk{index(1)+1,index(2)};
                nodec=blk{index(1)+2,index(2)};
                strb=['U_n' num2str(nodeb) '_0'];
                strc=['U_n' num2str(nodec) '_0'];
                % yout=char(yout,str, strb, strc);
                % y_type=[y_type 0 0 0];
                % nb_nodal_nodes=nb_nodal_nodes+3;
                % portType=[portType 1 1 1];
                if isempty(youtI)
                    youtI=char(str, strb, strc);
                else
                youtI=char(youtI,str, strb, strc);
                end
                y_typeI=[y_typeI 0 0 0];
                portTypeI=[portTypeI 1 1 1];
                nb_nodal_nodes=nb_nodal_nodes+3;
            else
                % if isempty(yout)
                %     yout=str;
                % else
                %     yout=char(yout,str);
                % end
                % y_type=[y_type 0];
                % nb_nodal_nodes=nb_nodal_nodes+1;
                % portType=[portType 1];

                if isempty(youtI)
                    youtI=str;
                else
                    youtI=char(youtI,str);
                end
                y_typeI=[y_typeI 0];
                portTypeI=[portTypeI 1];
                nb_nodal_nodes=nb_nodal_nodes+1;
            end
        end

        %%%%%%%%%%%%%%%55
        %%%%%attention
        % ces sources MATE ne sont pas encore incluses dans source!!!
        % ajouter ces sources ici.
        % ceci evite de prendre ces sorties dans yout!! (en les mettant
        % apres le calcul de yout)
        len=size(source,1);
        lenV=size(sourceV,1);
        lenI=size(sourceI,1);

        if MATEtype=='V'
            %from SPS doc: Voltage source: node1 is the positive terminal
            
            if m<=ind_sld
            %     source=[source; blk{index(1),4}   blk{index(1),3}   0 0 0 0];
            %     source=[source; blk{index(1)+1,4} blk{index(1)+1,3} 0 0 0 0];
            %     source=[source; blk{index(1)+2,4} blk{index(1)+2,3} 0 0 0 0];

                sourceV=[sourceV; blk{index(1),4}   blk{index(1),3}   0 0 0 0];
                sourceV=[sourceV; blk{index(1)+1,4} blk{index(1)+1,3} 0 0 0 0];
                sourceV=[sourceV; blk{index(1)+2,4} blk{index(1)+2,3} 0 0 0 0];
            else
                % source=[source; blk{index(1),3}   blk{index(1),4}   0 0 0 0];   % inverse node order sLD vs. 1ph
                sourceV=[sourceV; blk{index(1),3}   blk{index(1),4}   0 0 0 0];   % inverse node order sLD vs. 1ph
            end

            % sourcenames{len+1,1}=['MATEV_' blk{index(1),1}];
            sourcenamesV{lenV+1,1}=['MATEV_' blk{index(1),1}];
            if m<=ind_sld
            %     sourcenames{len+2,1}=['MATEV_' blk{index(1)+1,1}];
            %     sourcenames{len+3,1}=['MATEV_' blk{index(1)+2,1}];
                sourcenamesV{lenV+2,1}=['MATEV_' blk{index(1)+1,1}];
                sourcenamesV{lenV+3,1}=['MATEV_' blk{index(1)+2,1}];
            end
        else
            %from sps doc: Current source: Positive current flowing from node1 to node2 inside the source.
            
            if m<=ind_sld
                % source=[source; blk{index(1),3}   blk{index(1),4}   1 0 0 0];
                % source=[source; blk{index(1)+1,3} blk{index(1)+1,4} 1 0 0 0];
                % source=[source; blk{index(1)+2,3} blk{index(1)+2,4} 1 0 0 0];
                sourceI=[sourceI; blk{index(1),3}   blk{index(1),4}   1 0 0 0];
                sourceI=[sourceI; blk{index(1)+1,3} blk{index(1)+1,4} 1 0 0 0];
                sourceI=[sourceI; blk{index(1)+2,3} blk{index(1)+2,4} 1 0 0 0];
            else
                % source=[source; blk{index(1),4}   blk{index(1),3}   1 0 0 0];  % inverse node order sLD vs. 1ph
                sourceI=[sourceI; blk{index(1),4}   blk{index(1),3}   1 0 0 0];  % inverse node order sLD vs. 1ph
            end
            % sourcenames{len+1,1}=['MATEI_' blk{index(1),1}];
            sourcenamesI{lenI+1,1}=['MATEI_' blk{index(1),1}];
            if m<=ind_sld
                % sourcenames{len+2,1}=['MATEI_' blk{index(1)+1,1}];
                % sourcenames{len+3,1}=['MATEI_' blk{index(1)+2,1}];
                sourcenamesI{lenI+2,1}=['MATEI_' blk{index(1)+1,1}];
                sourcenamesI{lenI+3,1}=['MATEI_' blk{index(1)+2,1}];
            end
    end

    end
        %reorder MATE-block to put 3ph-SLD-V 1ph-V  3ph_SLD-I 1ph-I
    if isempty(yout)
        if isempty(youtI) && isempty(youtV)
            %yout=yout!
        elseif isempty(youtI)
            yout=char(youtV);
        elseif isempty(youtV)
            yout=char(youtI);
        else
            yout=char(youtI,youtV);
        end
    else
        if isempty(youtI) && isempty(youtV)
            %yout=yout!
        elseif isempty(youtI)
            yout=char(yout,youtV);
        elseif isempty(youtV)
            yout=char(yout,youtI);
        else
            yout=char(yout,youtI,youtV);
        end
    end
    y_type=[y_type y_typeI y_typeV];
    portType=[portType portTypeI portTypeV];

    source=[source; sourceI; sourceV];
    len=size(sourcenames,1);
    for k=1:size(sourcenamesI,1)
        sourcenames{len+k,1}=sourcenamesI{k,1};
    end
    len=size(sourcenames,1);
        for k=1:size(sourcenamesV,1)
        sourcenames{len+k,1}=sourcenamesV{k,1};
    end
    %sourcenames={sourcenames; sourcenamesI; sourcenamesV};

    % so, we now have the correct order of 'all current source' floowed by 'all
% voltage source'
% we must now, re-order in each of them according to blkmate node order.
orderV=[];
orderI=[];
allorder=[];
order=[];
for i=1:size(mateblk,1)
    if mateblk{i,3}=='V'
        orderV=[orderV str2num(mateblk{i,2})];
    end
    if mateblk{i,3}=='I'
        orderI=[orderI str2num(mateblk{i,2})];
    end
    %allorder=[allorder str2num(mateblk{i,2})+mateblk{i,3}*1];   %'V'*1=86   'I'*1= 73
end
%allorder=[orderI orderV+1000];  % the V and I source are laready ordered by type
allorder=[orderI orderV]; 
disp('MODIF MakeMAteInterface src_order 22 sept');
[a,src_order]=sort(allorder);
src_order
end



if isequal(step,'commentMATE')
    % comment the MATE sources
    for i=1:size(mateblk,1)
        set_param([mateblk{i,1} '/MATE_V'],'Commented','on');
        set_param([mateblk{i,1} '/MATE_I'],'Commented','on');
    end
    for i=1:size(mateblk_1ph,1)
        set_param([mateblk_1ph{i,1} '/MATE_V'],'Commented','on');
        set_param([mateblk_1ph{i,1} '/MATE_I'],'Commented','on');
    end
end

if isequal(step,'step1')   % not used anymore!
    for i=1:size(mateblk,1)
        blkh=get(get_param(mateblk{i,1},'Handle'));
        outblk=blkh.PortConnectivity(2).DstBlock(1);
        mateblk{i,2}=outblk;
        flag=0;
        for j=1:size(blk,1)
            if isequal( mateblk{i,2},get_param(blk{j,1},'Handle') )
                flag=j;
                type=get_param(mateblk{i,2},'MaskType');
                matetype=-1;
                if isequal (type, 'RLC (Three-Phase)')
                    rlctype=get_param(blk{j,1}, 'component_structure');

                    if isequal(rlctype,'ee.enum.rlc.structure.R')
                        matetype=0;  %V-type
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.L')
                        matetype=0;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.C')
                        matetype=1;  %I-type
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRL')
                        matetype=0;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRC')
                        matetype=0;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesLC')
                        matetype=0;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.SeriesRLC')
                        matetype=0;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRL')
                        matetype=1;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRC')
                        matetype=1;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelLC')
                        matetype=1;
                    elseif isequal(rlctype, 'ee.enum.rlc.structure.ParallelRLC')
                        matetype=1;
                    else
                    end
                    mateblk{i,3}=matetype;
                    % enable the MATE source
                    if matetype==0
                        set_param([mateblk{i,1} '/MATE_V'],'Commented','off');
                        set_param([mateblk{i,1} '/MATE_I'],'Commented','on');
                    elseif matetype==1
                        set_param([mateblk{i,1} '/MATE_V'],'Commented','on');
                        set_param([mateblk{i,1} '/MATE_I'],'Commented','off');
                    else
                    end
                else
                    error('MATE block error 1: only RLC SLD supported for now for MATE connection block')
                end
                break;
            end
        end
        if flag==0
            error('MATE block error 2')
        end
    end
end


end