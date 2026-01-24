function blk = AddRLsrcNodes(blk,blkmask,SLD_block_index,NodeNumber,blkgnd)
flag=0;
for i=1:size(blk,1)
    
    if blk{i,2}==SLD_block_index{6}(1)
        NodeNumber=NodeNumber+10000;

        paramNames=[];
        paramNames{1}='vline_rms';
        paramNames{2}='shift';
        paramNames{3}='freq';
        paramNames{4}='impedance_option';
        paramNames{5}='R';
        paramNames{6}='L';
        vals=ObtainParameterValue(blk{i,1},paramNames);
        ss=size(blk,1);
        if vals{4}>0
            % end rl block
            for k=1:size(blk,2)
                blk{ss+1,k}=blk{i,k};
            end
            blk{ss+1,2}=7; % rlc SLD
            no_n=blk{i,4};
            no_p=blk{i,3};
            newno=no_p+NodeNumber;
            blk{ss+1,3}=no_n;
            blk{ss+1,4}=newno;
            blk{i,4}=newno;
            blk{i,3}=no_p;
            flag=flag+1;
        end
    end

end
    if flag>0.5 %&& 1==0
        x=size(blk,1);
        y=size(blkgnd,1);
        blk=[blk(1:y,:); blk(x-flag+1:x,:); blk(y+1:x-flag,:)]; %put new rlc on top of blk{} afger gnd blocks
        disp(['block: ' blk(x-flag+1:x,1) 'put on top of blk list so to circonvent a curious power_statespace error']);
    end

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
        catch
            error(['3: not able to obtain params of block' blk])
        end
end
end