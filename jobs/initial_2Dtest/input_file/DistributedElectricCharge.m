%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Distributed loads whose direction is constant throughout the simulation
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function phi_distributed         =  DistributedElectricCharge(str)

phi_distributed                  =  zeros(str.mesh.volume.phi.n_nodes,1);
n_boundary_edges                 =  size(str.mesh.surface.phi.boundary_edges,2);
n_gauss                          =  size(str.quadrature.surface.bilinear.Chi,1);
w0                               =  1e-3; % Value of the applied traction
n_nodes_elem_boundary            =  size(str.mesh.surface.phi.boundary_edges,1);

ymin                             =  min(str.mesh.volume.phi.nodes(2,:));
for iedge=1:n_boundary_edges
    edge_nodes                   =  str.mesh.surface.phi.boundary_edges(:,iedge);
    xelem_boundary               =  str.mesh.volume.phi.nodes(:,edge_nodes);
    xcenter                      =  sum(xelem_boundary,2)/n_nodes_elem_boundary;
    Telem                        =  zeros(n_nodes_elem_boundary,1);
    if (xcenter(2) - ymin)<1e-5
       %-------------------------------------------------------------------
       % Compute kinematics to get DX_chi
       %-------------------------------------------------------------------
       kinematics                =  KinematicsFunctionSurface(str.geometry.dim,xelem_boundary,xelem_boundary,[0;0],...
                                                                str.fem.surface.bilinear.phi.N,str.fem.surface.bilinear.phi.DN_chi);         
       iso_Jacobian              =  kinematics.DX_chi_Jacobian;        
       weight                    =  str.quadrature.surface.bilinear.W_v;
       for igauss=1:n_gauss
           Telem                 =  Telem + w0*str.fem.surface.bilinear.phi.N(:,igauss)*(iso_Jacobian(igauss)*weight(igauss));
       end
       %-------------------------------------------------------------------
       % Assembly of the force vector
       %-------------------------------------------------------------------
       phi_nodes                 =  edge_nodes';
       phidofs                   =  phi_nodes;
       phi_distributed(phidofs)  =  phi_distributed(phidofs) + Telem(:);
    end
end









