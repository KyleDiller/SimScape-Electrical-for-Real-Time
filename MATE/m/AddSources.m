function [srcParam, src_Selector]=AddSources(source_parameter,nb_switch,nb_nodes)

nbs=size(source_parameter,1);
for i=1:nbs
    srcParam(1,i)=source_parameter{i,2}(1);
    srcParam(2,i)=source_parameter{i,2}(2);
    srcParam(3,i)=source_parameter{i,2}(3);
end
if nbs==0
    srcParam(:,1)=[-1 ; -1 ; -1];
end

if nb_switch>0
    src_Selector=[ones(1,nb_switch) 2:2+nbs-1 ones(1,nb_nodes)];
else
    src_Selector=[2:2+nbs-1 ones(1,nb_nodes)];
end



