%**************************************************************************
%**************************************************************************
%  1.This function plots the electric potential distribution 
%**************************************************************************
%**************************************************************************

function    plot_phi(str)   
orderN     =  [1 2 3 4];   
%------------------------------------------
% Plot electric potential
%------------------------------------------
figure(2000)
for ielem=1:str.n_elem
        order(1:length(orderN))              =   str.connectivity(ielem,orderN);
        phielem                              =   str.phi(order);                
       xplot                                =   str.Eulerian_x(:,order);        
%        xplot                                =   str.Lagrangian_X(:,order);        
%        hv=fill(xplot(1,:),xplot(2,:),phielem,'EdgeColor','none')  ;
        hv=fill(xplot(1,:),xplot(2,:),phielem)  ;
    hold on
end  
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;
end
