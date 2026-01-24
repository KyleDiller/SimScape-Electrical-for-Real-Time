function [blk, Vnode,cur_br_all, sensorV,sensorI] = MergeCurrentSensorNodes(blk,SLD_block_index)
Vnode=[];

cur_br_all={};
sensorV={};
vv=0;
sensorI={};
ii=0;


k=1;
for i=1:size(blk,1)
    if blk{i,2}==SLD_block_index{4}(1) || blk{i,2}==SLD_block_index{4}(2)
        cur_br_right=[];
        cur_br_left=[];
        signl=[];
        signr=[];
        if blk{i,2}==SLD_block_index{4}(1)
            ii=ii+1;
            sensorI{ii,1}=blk{i,1};
        elseif blk{i,2}==SLD_block_index{4}(2)  % this is a 3 phase sensor
     
            if isequal(blk{i,1}(end-1:end),'_B')
                iii=i-1;
            elseif isequal(blk{i,1}(end-1:end),'_C')
                iii=i-2;
            else
                iii=i;
            end
            ii=ii+1;
            sensorI{ii,1}=blk{iii,1};
            vv=vv+1;
            sensorV{vv,1}=blk{iii,1};        
            
        else
            error('Unknown sensor...')
        end

        staynode=blk{i,4};
        searchnode=blk{i,3};
        if searchnode==0
            searchnode=staynode;
            staynode=0;
            disp('Current snesor connectet to ground: PLEASE TEST!')
        end

        for j=1:size(blk,1)
            if isequal(blk{j,3},searchnode)
                disp(['1 search node: ' num2str(searchnode) ' staynode:' num2str(staynode)]);
                if isempty(find(abs(blk{j,2})==SLD_block_index{3})) && j~=i % forget phase splitter blocks Note that there type
                    % have been negated in the reduction
                    % routine that is why the abs()
                    if isempty(find(abs(blk{j,2})==SLD_block_index{5}))  % is rlc or source?
                        cur_br_left=[cur_br_left j];
                        signl=[signl '+'];    % +-  scr
                    else
                        cur_br_right=[cur_br_right j];
                        signr=[signr '-'];  %rlc
                    end
                    blk{j,3}=staynode;
                end
            end
            if isequal(blk{j,4},searchnode)
                if isempty(find(abs(blk{j,2})==SLD_block_index{3})) && j~=i
                    disp(['1 search node: ' num2str(searchnode) ' staynode:' num2str(staynode)]);
                    if isempty(find(abs(blk{j,2})==SLD_block_index{5}))  % is rlc or source?
                        cur_br_left=[cur_br_left j];
                        signl=[signl '-'];
                    else
                        cur_br_right=[cur_br_right j];
                        signr=[signr '+'];
                    end
                    blk{j,4}=staynode;
                end
            end
        end
        if blk{i,2}==SLD_block_index{4}(2)
            Vnode=[Vnode staynode];
        end
        ind=find(Vnode==searchnode);
        if ~isempty(ind)  % this means that a new current
            % sensors is merging nodes and the disaprearing node
            % has already been registered in Vnode
            % so we need to update it
            Vnode(ind)=staynode;
        end
        cur_br_all{k,1}=blk{i,1};
        cur_br_all{k,2}=cur_br_right;   %rlc
        cur_br_all{k,3}=cur_br_left;    %source
        cur_br_all{k,4}=signr;
        cur_br_all{k,5}=signl;
        k=k+1;

    end


end

end  % end function