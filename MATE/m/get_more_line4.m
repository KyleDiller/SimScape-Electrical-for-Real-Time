function [new_line_handle, new_block_handle, new_port_handle,ssc_flag] = get_more_line4(BlockHandle, PortHandle,silent)
%UNTITLED Summary of this function goes here

new_line_handle=[];
new_block_handle=[];
new_port_handle=[];
ssc_flag=0;

typ=get_param(BlockHandle,'BlockType');

if isequal(typ,'SimscapeBlock')
    param=get(BlockHandle);
    Allports=[param.PortHandles.LConn param.PortHandles.RConn];
    portnumber=find(PortHandle==Allports);
    if silent
    disp(['4.1: Simscape block found ' getfullname(BlockHandle) ' Port: ' num2str(portnumber)]); 
    end
    new_block_handle=[BlockHandle];
    new_port_handle=[PortHandle];
    ssc_flag=1;
    return
elseif isequal(typ,'PMIOPort')

    portnumber=str2num(get_param(BlockHandle,'Port')); % port number is easy here
   
    par=get_param(BlockHandle,'Parent');
    blks=get_param(par,'Blocks');
    parh=get_param(par,'Handle');
    pc=get_param(par,'PortConnectivity');
    ThePort=pc(portnumber);  % upper subsystem level port
    pr=get(parh);
    
        % find the index of the port in Lconn Rconn combined list   
   %  allportshandle=[pr.PortHandles.LConn pr.PortHandles.RConn];
   %  upperblk=allportshandle(portnumber);
   %  upbh=get(upperblk);
   % length(upbh.PortConnectivity.DstPort)
   for k=1:length(ThePort.DstBlock)
       if isequal( get_param( getfullname(ThePort.DstBlock(k)) ,'Commented'),'off')
           new_block_handle=[new_block_handle ThePort.DstBlock(k)];
           new_port_handle=[new_port_handle ThePort.DstPort(k)];
           if silent
           disp(['PMIO ' getfullname(BlockHandle) ' port:' num2str(portnumber) '->' getfullname(ThePort.DstBlock(k))  'port: ' num2str(k)]);
           end
       end
   end
   allportshandle=[pr.PortHandles.LConn pr.PortHandles.RConn];
   new_line_handle=[allportshandle(portnumber)];
   
 
  
elseif isequal(typ,'SubSystem')
     % find the index of the port in Lconn Rconn combined list
     pr=get(BlockHandle);
     allportshandle=[pr.PortHandles.LConn pr.PortHandles.RConn];
     for k=1:length(allportshandle)
         if allportshandle(k)==PortHandle
             break;  % k is the indice of the subsystem port (seen from oput of the subsystem)
         end
     end
     subblks=get_param(BlockHandle,'Blocks');
     blkname=getfullname(BlockHandle);

     for ii=1:size(subblks,1)
         if isequal(get_param([blkname '/' subblks{ii}],'BlockType'),'PMIOPort')
                portnumber=str2num(get_param([blkname '/' subblks{ii}],'Port'));
                if k==portnumber
                    % on a trouve le PMIOport dans le subsystem!
                    % c'est le ii'th block
                    break;
                end
         end
     end
     tmp=get_param([blkname '/' subblks{ii}], "LineHandles");
     new_line_handle=[tmp.LConn tmp.RConn];
     %get_connectivity here
     blkh=get(get_param([blkname '/' subblks{ii}],"Handle"));
     for j=1: length(blkh.PortConnectivity.DstPort)
         if isequal(get_param([blkname],'Commented'),'off')
             new_block_handle=[new_block_handle blkh.PortConnectivity.DstBlock(j)];
             new_port_handle=[new_port_handle blkh.PortConnectivity.DstPort(j)];
         end
     end
     if silent
     disp(['Subsystem found 6.1 connection: ' [blkname '/' subblks{ii}] ' Line handle: ' num2str(new_line_handle)]);
     end
end

end