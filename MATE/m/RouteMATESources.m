function RouteMATESources(modele,order,iem)
%(modele, src_order,'InputElectricMixer');
%to order the MATE source correctly.
%blk='RouteTestModel_11/X'
%order=[3     1     2     4     5     6]
%phase=6;
% 'a'=97
% 'A'=65
par=get_param(modele,'Parent');
mod=[par '/' iem];


phase=length(order);

lin = find_system(mod, 'FindAll', 'on', 'type', 'line');
delete_line(lin)

for i=1:phase
    h1=get_param([mod, '/' [char(i+64) ]],'PortHandles');
    h2=get_param([mod, '/' [char(order(i)+96) ]],'PortHandles');
    add_line(mod,h1.RConn(1),h2.RConn(1));
end


