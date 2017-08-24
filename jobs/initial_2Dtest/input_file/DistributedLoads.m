%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Distributed loads whose direction is constant throughout the simulation
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function P_distributed       =  DistributedLoads(str)

P_distributed                =  zeros(str.geometry.dim*str.mesh.volume.x.n_nodes,1);
n_boundary_edges             =  size(str.mesh.surface.x.boundary_edges,2);
n_gauss                      =  size(str.quadrature.surface.bilinear.Chi,1);
p0                           =  zeros(3,1); % Value of the applied traction
n_nodes_elem_boundary        =  size(str.mesh.surface.x.boundary_edges,1);

ymin                         =  min(str.mesh.volume.x.nodes(2,:));
for iedge=1:n_boundary_edges
    edge_nodes               =  str.mesh.surface.x.boundary_edges(:,iedge);
    xelem_boundary           =  str.mesh.volume.x.nodes(:,edge_nodes);
    xcenter                  =  sum(xelem_boundary,2)/n_nodes_elem_boundary;
    Telem                    =  zeros(str.geometry.dim,n_nodes_elem_boundary);
    if (xcenter(2) - ymin)<1e-5
       %-------------------------------------------------------------------
       % Compute kinematics to get DX_chi
       %-------------------------------------------------------------------
       kinematics            =  KinematicsFunctionSurface(str.geometry.dim,xelem_boundary,xelem_boundary,[0;0;0],...
                                                          str.fem.surface.bilinear.x.N,str.fem.surface.bilinear.x.DN_chi);         
       iso_Jacobian          =  kinematics.DX_chi_Jacobian;        
       weight                =  str.quadrature.surface.bilinear.W_v;
       for igauss=1:n_gauss
           Telem             =  Telem + (p0*str.fem.surface.bilinear.x.N(:,igauss)')*(iso_Jacobian(igauss)*weight(igauss));
       end
       %-------------------------------------------------------------------
       % Assembly of the force vector
       %-------------------------------------------------------------------
       x_nodes               =  edge_nodes';
       xdofs                 =  zeros(3,size(x_nodes,2));
       for idim=1:3
           xdofs(idim,:)     =  (x_nodes-1)*3 + idim;
       end
       xdofs                 =  reshape(xdofs,size(xdofs,1)*size(xdofs,2),1);
       P_distributed(xdofs)  =  P_distributed(xdofs) + Telem(:);
    end
end




