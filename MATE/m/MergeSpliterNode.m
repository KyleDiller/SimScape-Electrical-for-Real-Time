function [blk] = MergeSpliterNode(blk,blkmask,SLD_block_index, offsetnode)

for i=1:size(blk,1)
    if blk{i,2}==SLD_block_index{3}
        if isempty(blk{i,5})  % to find the single-phase links that replace the 3-phase splitter
            staynode=blk{i,4};
            searchnode=blk{i,3};
            for j=1:size(blk,1)
                if isequal(blk{j,3},searchnode)
                    blk{j,3}=staynode;
                end
                if isequal(blk{j,4},searchnode)
                    blk{j,4}=staynode;
                end
            end
            blk{i,2}=0-blk{i,2};   % negate the branch to mark it
        end
    end
end


end