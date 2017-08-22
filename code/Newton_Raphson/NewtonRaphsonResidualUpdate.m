%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function updates the forces and matrices for the next Newton Raphson 
% iteration.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str              =  NewtonRaphsonResidualUpdate(str)

%--------------------------------------------------------------------------
% Compute Neumann forces/electric charges
%--------------------------------------------------------------------------
str                       =  NeumannBcs(str);
%--------------------------------------------------------------------------
% Initialise the residual
%--------------------------------------------------------------------------
str.assembly.Residual     =  str.assembly.total_force;
dofs                      =  str.mesh.dimesion.volume.x.n_nodes + str.mesh.volume.phi.n_nodes;        
%--------------------------------------------------------------------------
% Compute the residual for static or dynamic simulations
%--------------------------------------------------------------------------
switch str.data.analysis
    case 'static'
        str.assembly.Residual(dofs,...
            1)            =  str.assembly.total_force(dofs) -  ...
                             str.bc.Neumann.force_vector*str.NR.accumulated_factor;
    case 'dynamic'
        switch str.time_integrator.type
            case 'generalised_alpha'
                 str.assembly.Residual(dofs,...
                     1)   =  str.assembly.total_force(dofs,1) -...
                             ((1-time_integrator.alpha)*str.bc.Neumann.force_vector + ...
                             str.time_integrator.alpha*str.bc.Neumann.force_vector);
            case 'Newmark_beta'
                str.assembly.Residual(dofs,...
                    1)    =  str.assembly.total_force(dofs) -  str.bc.Neumann.force_vector;                
        end
end


