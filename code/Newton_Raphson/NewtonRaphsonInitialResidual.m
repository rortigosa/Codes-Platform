%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Compute initial residual taking into account Dirichlet boundary
% conditions
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                           =  NewtonRaphsonInitialResidual(str,old_solution)

%--------------------------------------------------------------------------
% Compute Neumann forces/electric charges
%--------------------------------------------------------------------------
str                                    =  NeumannBcs(str);
%--------------------------------------------------------------------------
% Initial residual for static or dynamic cases
%--------------------------------------------------------------------------
switch str.data.analysis
    case 'static'
         %-----------------------------------------------------------------
         % Difference between the Neumann forces vectors between each load increment
         %-----------------------------------------------------------------
         former_force_vector           =  (str.NR.accumulated_factor - str.NR.load_factor)*str.bc.Neumann.force_vector;
         str.bc.Neumann.force_vector   =  str.NR.accumulated_factor*str.bc.Neumann.force_vector;
         diff_force_vector             =  str.bc.Neumann.force_vector - former_force_vector;
         %-----------------------------------------------------------------
         % Compute the residual based on the Neumann forces
         %-----------------------------------------------------------------        
         dofs                          =  str.geometry.dim*str.mesh.volume.x.n_nodes + str.mesh.volume.phi.n_nodes;
         str.assembly.Residual(1:dofs,...
             1)                        =  str.assembly.Residual(1:dofs,1) - diff_force_vector;
         diff_fields                   =  DiffFieldsFunction(str.data.formulation,str.solution,old_solution);
         %-----------------------------------------------------------------
         % Stiffness matrix K(free,constrained).
         %-----------------------------------------------------------------
         constrained_dofs              =  (1:size(diff_fields,1))';
         constrained_dofs(str.bc.Dirichlet.freedof,...
             1)                        =  [];
         Delta_constrained_dofs        =  diff_fields(constrained_dofs);
         reduced_matrix                =  str.assembly.K_total(str.bc.Dirichlet.freedof,constrained_dof);
         %-----------------------------------------------------------------
         % Residual=.K(free,constrained)*Deltaconstrained.
         %-----------------------------------------------------------------
         str.assembly.Residual(str.bc.Dirichlet.freedof,...
             1)                        =  str.assembly.Residual(str.bc.Dirichlet.freedof,1) + ...
                                          reduced_matrix*Delta_constrained_dofs;
end
end