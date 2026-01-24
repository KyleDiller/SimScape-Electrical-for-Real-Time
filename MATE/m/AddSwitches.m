function [swParam,masterSwSelOrder]=AddSwitches(modele,switchnames, swmixerblk)

msi=find_system(modele,'FindAll','on','MaskType','MATE Switch Identifier')

for i=1:size(msi,1)
    fullnam=getfullname(msi(i));
    tag=get_param([fullnam '/Goto'],'Gototag');
    swn=get_param(msi(i),'bn');
    flag=0;
    for j=1:size(switchnames,1)
        swn2=get_param(switchnames{j,1},'Name');
        if isequal(swn2, swn)
            if flag==0
            indice=j;
            flag=flag+1;
            end
        end
        if isequal(swn2, swn)
            flag=flag+1;    % check for _B and _C  (A B C always follow each other)
        end
    end
    if flag==0
        error('AddSwitches code 1')
    end
    swParam{i,1}=tag;
    swParam{i,2}=indice;
    swParam{i,3}=flag-1; % number of expended gate signals
end

par=get_param(modele,'Parent');
masterSwSel=[];
for i=1:size(swParam,1)
    set_param([par '/' swmixerblk '/From' num2str(i)], 'Gototag',swParam{i,1})
    masterSwSel=[masterSwSel  swParam{i,2}];
    str='[';
    for j=1:swParam{i,3}
        str=[str ' 1 '];
    end
    str=[str ']' ];
    set_param([par '/' swmixerblk '/Selector' num2str(i)],'Indices', str)
end
[tmp,masterSwSelOrder]=sort(masterSwSel);
for i=1:size(swParam,1)
    set_param([par '/' swmixerblk '/MainSelector' num2str(i)],'Indices', ['[' num2str(masterSwSelOrder(i)) ']']);
end


% nbs=size(source_parameter,1);
% for i=1:nbs
%     srcParam(1,i)=source_parameter{i,2}(1);
%     srcParam(2,i)=source_parameter{i,2}(2);
%     srcParam(3,i)=source_parameter{i,2}(3);
% end
% if nb_switch>0
%     src_Selector=[ones(1,nb_switch) 2:2+nbs-1 ones(1,nb_nodes)]
% else
%     src_Selector=[2:2+nbs-1 ones(1,nb_nodes)]
% end



