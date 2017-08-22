%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Add to the stiffness matrix contributions from inertia, damping, internal
% work...
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function str                                 =  StiffnessMatrixTotalAssembly(str)

mec_dofs                                     =  size(str.solution.x.Eulerian_x(:),1);
%--------------------------------------------------------------------------
% Multiply matrices by the time integration factors.
%--------------------------------------------------------------------------
switch str.data.analysis
    case 'static'
        %------------------------------------------------------------------
        % Final stiffness matrix for statics.
        %------------------------------------------------------------------
        str.assembly.total_matrix            =  str.assembly.K_total;
    case 'dynamic'
        switch str.time_integrator.type
            case 'Newmark_beta'
                 str.assembly.total_matrix   =  str.assembly.K_total;
                 Dt                          =  str.time_integrator.Deltat;
                 beta                        =  str.time_integrator.beta;
                 gamma                       =  str.time_integrator.gamma;
                 
                 Mfactor                     =  1/(beta*Dt^2);
                 Dfactor                     =  gamma/(beta*Dt);
                 str.assembly.total_matrix(1:mec_dofs,...
                     1:mec_dofs)             =  Mfactor*str.assembly.M_total + ...
                                                Dfactor*str.assembly.Damping_total + ...
                                                str.assembly.K_total(1:mec_dofs,1:mec_dofs);
        end
end
