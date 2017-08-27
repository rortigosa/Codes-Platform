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
p0                           =  zeros(str.geometry.dim,1); % Value of the applied traction
n_nodes_elem_boundary        =  size(str.mesh.surface.x.boundary_edges,1);

xmax                         =  max(str.mesh.volume.x.nodes(1,:));
for iedge=1:n_boundary_edges
    edge_nodes               =  str.mesh.surface.x.boundary_edges(:,iedge);
    xelem_boundary           =  str.mesh.volume.x.nodes(:,edge_nodes);
    xcenter                  =  sum(xelem_boundary,2)/n_nodes_elem_boundary;
    Telem                    =  zeros(str.geometry.dim,n_nodes_elem_boundary);
    if abs(xcenter(1) - xmax)<1e-5
       %-------------------------------------------------------------------
       % Compute kinematics to get DX_chi
       %-------------------------------------------------------------------
       kinematics            =  KinematicsFunctionSurface(str.geometry.dim,xelem_boundary,xelem_boundary,zeros(str.geometry.dim,1),...
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
       xdofs                 =  zeros(str.geometry.dim,size(x_nodes,2));
       for idim=1:str.geometry.dim
           xdofs(idim,:)     =  (x_nodes-1)*str.geometry.dim + idim;
       end
       xdofs                 =  reshape(xdofs,size(xdofs,1)*size(xdofs,2),1);
       P_distributed(xdofs)  =  P_distributed(xdofs) + Telem(:);
    end
end




