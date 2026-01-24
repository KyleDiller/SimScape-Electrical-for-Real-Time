function [outstr,y_type,output_sig] = BuildSPSYout(Vnode, switches,cur_br_all,rlcnames,sourcenames,switchnames)
% % The measurements that you request follow, in any order.
% %
% y_u1='U_n10_11';			%U_Sw1= Voltage across Sw1
% y_u2='U_n11_12';			%U_Sw2= Voltage across Sw2
% y_i3='I1';			%I1= Switch current Sw1
% y_i4='I2';			%I2= Switch current Sw2
% y_u5='U_n2.1_0';			%U_sat= Voltage across saturable reactor 
% y_i6='I_b1';			%I1 measurement
% y_u7='U_n11_0';			%V2 measurement
% y_u8='U_n12_0';			%V3 measurement
% 
% yout=char(y_u1,y_u2,y_i3,y_i4,y_u5,y_i6,y_u7,y_u8);								% outputs
% % y_type=[0,0,1,1,0,1,0,0];				%output types; 0=voltage 1=current

outstr=[];
y_type=[];
output_sig={};
for i=1:size(switches,1)
    str='U_n';
    str=[str num2str(switches(i,1)) '_' num2str(switches(i,2))];
    if isempty(outstr)
        outstr=str;
    else
        outstr=char(outstr, str);
    end
    y_type=[y_type 0];
    output_sig{i}=switchnames{i};
end
off=size(switches,1)
for i=1:length(Vnode)
    str='U_n';
    str=[str num2str(Vnode(i)) '_0'];
    if isempty(outstr)
        outstr=str;
    else
        outstr=char(outstr, str);
    end
    y_type=[y_type 0];
    output_sig{i+off}=str;
end

for i=1:size(cur_br_all,1)
    % get branches blk index
    str=[];
    for j=1:length(cur_br_all{i,2})

        %cur_br_all{i,2}(j); % block number for rlc
        flag=0;
        for k=1:size(rlcnames,1)
            if rlcnames{k,2}==cur_br_all{i,2}(j)
                %str=[str num2str(k)];
                flag=1;
                break;
            end
        end
        if flag==0
            error('error BuildSPSYout 1')
        end
        if j==1
            str=[cur_br_all{i,4}(j) 'I_b' num2str(k)];
        end
        if j>1
            str=[str cur_br_all{i,4}(j) 'I_b' num2str(k)];  % cur_br_all{i,4}(j) is the sign
        end
        
    end
    for j=1:length(cur_br_all{i,3})


        %cur_br_all{i,3}(j); % block number for source
        flag=0;
        for k=1:size(sourcenames,1)
            if sourcenames{k,2}==cur_br_all{i,3}(j)
                %str=[str num2str(k)];
                flag=1;
                break;
            end
        end
        if flag==0
            error('error BuildSPSYout 2')
        end
        if j==1 && isempty(str)
            %str='I'
            str=[cur_br_all{i,5}(j) 'I' num2str(k)];
        else
            %str=[str '+I'];
            str=[str cur_br_all{i,4}(j) 'I' num2str(k)];
        end
        
    end
    if ~isempty(str)
    outstr=char(outstr, str);
    y_type=[y_type 1];
    end
   
end

end