function [newblk] = ExpendSLDNodes(blk,blkmask,SLD_block_index, offsetnode)


newblk=cell(0,size(blk,2));
% expand std SLD blocks

for i=1:size(blk,1)
    if isequal( blkmask{blk{i,2},4},-1)
        newblk=[newblk ; cell(1,size(blk,2))];
        for j=1:size(blk,2)
            newblk{end,j}= blk{i,j};
        end
    else
        if ~isempty( find( SLD_block_index{1}==blk{i,2} ) )  % this is an SLD format
            newblk=[newblk ; cell(1,size(blk,2))];
            if ~isempty( find( SLD_block_index{6}==blk{i,2} ) )  % SLD source with RL (or not!)
                
            end
            for j=1:size(blk,2)
                newblk{end,j}= blk{i,j}; % copy phase A as is
            end
            for j=1:length(SLD_block_index{1})
                if blk{i,2}==SLD_block_index{1}(j)
                    newblk=[newblk ; cell(1,size(blk,2))];  % phase B following with offssetnode node number
                    newblk{end,1}=[blk{i,1} '_B'];
                    newblk{end,2}=blk{i,2};
                    if blk{i,3}<0.5
                        newblk{end,3}=blk{i,3};
                    else
                        newblk{end,3}=blk{i,3}+offsetnode;
                    end
                    if blk{i,4}<0.5
                        newblk{end,4}=blk{i,4};
                    else
                        newblk{end,4}=blk{i,4}+offsetnode;
                    end
                    newblk=[newblk ; cell(1,size(blk,2))];  % phase B following with offssetnode node number
                    newblk{end,1}=[blk{i,1} '_C'];
                    newblk{end,2}=blk{i,2};
                    if blk{i,3}<0.5
                        newblk{end,3}=blk{i,3};
                    else
                        newblk{end,3}=blk{i,3}+2*offsetnode;
                    end
                    if blk{i,4}<0.5
                        newblk{end,4}=blk{i,4}; % case of 0 ground.
                    else
                        newblk{end,4}=blk{i,4}+2*offsetnode;
                    end
                    %break; % next blk
                end
            end
        end
        %disp(['blk: ' blk{i,1} ' type: ' blk{i,2}])
        if ~isempty( find( SLD_block_index{2}==blk{i,2} ) )   % this is an 3ph format
            %disp(['HIT@@@@ blk: ' blk{i,1} ' type: ' blk{i,2}])
            newblk=[newblk ; cell(1,size(blk,2))];
            newblk{end,1}=blk{i,1};  %leave A name unchanged
            newblk{end,2}=blk{i,2};
            newblk{end,3}=blk{i,3};
            newblk{end,4}=blk{i,6};
            newblk=[newblk ; cell(1,size(blk,2))];
            newblk{end,1}=[blk{i,1} '_B'];
            newblk{end,2}=blk{i,2};
            newblk{end,3}=blk{i,4};
            newblk{end,4}=blk{i,7};
            newblk=[newblk ; cell(1,size(blk,2))];
            newblk{end,1}=[blk{i,1} '_C'];
            newblk{end,2}=blk{i,2};
            newblk{end,3}=blk{i,5};
            newblk{end,4}=blk{i,8};
            %break; %next block
        end
    end
end



% expend phase-splitters
% first LConn is to be expended only
% the 3-phase side has already been numbered
% so we create 3 single phase links
% continue with previous k index

    for i=1:size(blk,1)
        if blk{i,2}==SLD_block_index{3}

            newblk=[newblk ; cell(1,size(blk,2))]; 
            newblk{end,1}=[blk{i,1} '_A'];
            newblk{end,2}=blk{i,2};

            newblk{end,3}=blk{i,3};

            newblk{end,4}=blk{i,4};   % node already attributed on the 3ph side

            newblk=[newblk ; cell(1,size(blk,2))]; 
            newblk{end,1}=[blk{i,1} '_B'];
            newblk{end,2}=blk{i,2};
            if blk{i,3}<0.5
                newblk{end,3}=blk{i,3};
            else
                newblk{end,3}=blk{i,3}+offsetnode;
            end

            newblk{end,4}=blk{i,5}; % node already attributed on the 3ph side

            newblk=[newblk ; cell(1,size(blk,2))]; 
            newblk{end,1}=[blk{i,1} '_C'];
            newblk{end,2}=blk{i,2};
            if blk{i,3}<0.5
                newblk{end,3}=blk{i,3};
            else
                newblk{end,3}=blk{i,3}+2*offsetnode;
            end

            newblk{end,4}=blk{i,6}; % node already attributed on the 3ph side
        end
    end

end