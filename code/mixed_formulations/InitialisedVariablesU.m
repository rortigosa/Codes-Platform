function solution                   =  InitialisedVariablesU(geom,mesh)

%--------------------------------------------------------------------------
%  Kinematics
%--------------------------------------------------------------------------
solution.x.Lagrangian_X             =  mesh.volume.x.nodes;
solution.x.Eulerian_x               =  solution.x.Lagrangian_X;
solution.x.Eulerian_x_dot           =  zeros(geom.dim,mesh.volume.x.n_nodes);
solution.x.Eulerian_x_dot_dot       =  zeros(geom.dim,mesh.volume.x.n_nodes);
solution.x.Eulerian_x_dot_dot_old   =  zeros(geom.dim,mesh.volume.x.n_nodes);
solution.x.Eulerian_x_old           =  solution.x.Lagrangian_X;
%--------------------------------------------------------------------------
% Total number of dofs
%--------------------------------------------------------------------------
solution.n_dofs                     =  geom.dim*mesh.volume.x.n_nodes;

                                   
                                   