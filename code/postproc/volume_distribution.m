function [str] =  volume_distribution(str)

str.volume                     =  zeros(str.n_elem,1);
str.volume0                    =  zeros(str.n_elem,1);
for ielem=1:str.n_elem
    nodes_elem                 =  str.connectivity(ielem,:);
    xelem                      =  str.Eulerian_x(:,nodes_elem);    
    Xelem                      =  str.Lagrangian_X(:,nodes_elem);
    phielem                    =  zeros(str.n_node_elem,1);
    [str]                      =  gradients(xelem,Xelem,phielem,str);     
    for igauss=1:size(str.quadrature.Chi,1)
        J                      =  str.grad.J(igauss,1);
        DX_chi                 =  det(str.grad.DX_chi(:,:,igauss));        
        W                      =  str.quadrature.W_v(igauss);
        str.volume(ielem,1)    =  str.volume(ielem,1)  +  J*DX_chi*W;
        str.volume0(ielem,1)   =  str.volume0(ielem,1) +  DX_chi*W;
    end        
end

str.volume_ratio               =  str.volume./str.volume0;

%------------------------------------------------------
%  Plot specified distribution.
%------------------------------------------------------
orderN     =  [1 2 3 4];   
figure(3000)
for ielem=1:str.n_elem
        order(1:length(orderN))              =   str.connectivity(ielem,orderN);
        phielem                              =   str.volume_ratio(ielem);                
        xplot                                =   str.Eulerian_x(:,order);        
%        xplot                               =   str.Lagrangian_X(:,order);        
        hv=fill(xplot(1,:),xplot(2,:),phielem,'EdgeColor','none')  ;
 %       hv=fill(xplot(1,:),xplot(2,:),phielem)  ;
    hold on
end  
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;
