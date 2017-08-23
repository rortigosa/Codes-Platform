function element_assembly    =  ElementResidualInitialisationFormulationBoundaryBEM(dim,mesh)

ndof_x                       =  dim*mesh.volume.x.n_node_elem;
ndof_phi                     =  mesh.volume.phi.n_node_elem;
ndof_q                       =  mesh.surface.q.n_node_elem;
%--------------------------------------------------------------------------
% Residual vectors per element
%--------------------------------------------------------------------------
element_assembly.Tq          =  zeros(ndof_q,1);  
%--------------------------------------------------------------------------
% Stiffness matrices per element
%--------------------------------------------------------------------------
element_assembly.Kqx         =  zeros(ndof_q,ndof_x);
element_assembly.Kqphi       =  zeros(ndof_q,ndof_phi);
element_assembly.Kqq         =  zeros(ndof_q,ndof_q);
element_assembly.Kqphiprime  =  zeros(ndof_q,ndof_phi);
