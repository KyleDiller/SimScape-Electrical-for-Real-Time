function  netlistprint(blk)
%printout
nb_blk=size(blk,1);
for iii=1:nb_blk
    if iii==23
        ssss=1;
    end
    str=[];
    ind=find(blk{iii,1}=='/');
    nam=blk{iii,1}(ind(end)+1:end);
    str=[str nam  ' Node ' num2str(blk{iii,3}) ' node ' num2str(blk{iii,4})];
    kkk=5;
    while  size(blk,2)>=kkk
        if (blk{iii,kkk}>-0.5)
            str=[str ' Node ' num2str(blk{iii,kkk})];
            if kkk>size(blk,2)
                break;
            end
        end
        kkk=kkk+1;
    end
    disp(str);
end
end