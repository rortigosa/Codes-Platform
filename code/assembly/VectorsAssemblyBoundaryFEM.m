function assembly                   =  VectorsAssemblyBoundaryFEM(iedge,dim,mesh,element_assembly,assembly)

%-------------------------------------------------------------------------- 
% x
%--------------------------------------------------------------------------
n_nodes_elem                        =  size(mesh.surface.x.boundary_edges,1);
T1                                  =  repmat(1:dim,n_nodes_elem,1);
T1                                  =  T1';
T1                                  =  reshape(T1,dim*n_nodes_elem,1);

nodes                               =  mesh.surface.x.boundary_edges(:,iedge);
T2                                  =  ((nodes - 1)*dim);
T2                                  =  repmat(T2,1,dim);
T2                                  =  T2';
T2                                  =  reshape(T2,dim*n_nodes_elem,1);

global_x_dof                        =  T1 + T2;
assembly.Tx(global_x_dof,1)         =  assembly.Tx(global_x_dof,1) + element_assembly.Tx;
%--------------------------------------------------------------------------
%  phi
%--------------------------------------------------------------------------
global_phi_dof                      =  mesh.surface.phi.boundary_edges(:,iedge);
assembly.Tphi(global_phi_dof,1)     =  assembly.Tphi(global_phi_dof,1) + element_assembly.Tphi;         
%--------------------------------------------------------------------------
%  Pressure
%--------------------------------------------------------------------------
global_q_dof                       =  mesh.surface.q.connectivity(:,iedge);
assembly.Tq(global_q_dof,1)        =  assembly.Tq(global_q_dof,1) + element_assembly.Tq;         
