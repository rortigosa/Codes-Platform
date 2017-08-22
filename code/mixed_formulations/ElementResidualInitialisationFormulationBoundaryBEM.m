function str                              =  ElementResidualInitialisationFormulationBoundaryBEM(str)

ndof_x                                    =  mesh.dim*mesh.volume.x.n_node_elem;
ndof_phi                                  =  mesh.volume.phi.n_node_elem;
ndof_q                                    =  mesh.surface.q.n_node_elem;
%--------------------------------------------------------------------------
% Residual vectors per element
%--------------------------------------------------------------------------
str.assembly.element_assembly.Tq          =  zeros(ndof_q,1);  
%--------------------------------------------------------------------------
% Stiffness matrices per element
%--------------------------------------------------------------------------
str.assembly.element_assembly.Kqx         =  zeros(ndof_q,ndof_x);
str.assembly.element_assembly.Kqphi       =  zeros(ndof_q,ndof_phi);
str.assembly.element_assembly.Kqq         =  zeros(ndof_q,ndof_q);
str.assembly.element_assembly.Kqphiprime  =  zeros(ndof_q,ndof_phi);
