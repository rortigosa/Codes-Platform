function solution         =  InitialisedVariablesElectroBEMOnlyElectro(geom,mesh,quadrature)

%--------------------------------------------------------------------------
% Kinematics
%--------------------------------------------------------------------------
solution.x.Lagrangian_X   =  mesh.volume.x.nodes;
%--------------------------------------------------------------------------
% Electric potential
%--------------------------------------------------------------------------
solution.phi              =  zeros(mesh.volume.phi.n_nodes,1);
%--------------------------------------------------------------------------
% Electric flux in the boundary
%--------------------------------------------------------------------------
solution.q                          =  zeros(mesh.surface.q.n_nodes,1);
%--------------------------------------------------------------------------
% Electric displacement field
%--------------------------------------------------------------------------
solution.D0               =  zeros(geom.dim,size(quadrature.volume.bilinear.Chi,1),mesh.volume.n_elem);
%--------------------------------------------------------------------------
% Total number of dofs
%--------------------------------------------------------------------------
solution.n_dofs                     =  mesh.volume.phi.n_nodes + ...
                                       mesh.surface.q.n_nodes;
