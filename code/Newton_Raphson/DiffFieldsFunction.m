%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function computes the difference between the fields of a specific
% formulation at two different stages of the simulation
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function diff_fields  =  DiffFieldsFunction(formulation,solution,old_solution)

switch formulation
    case 'electro'
        diff_x        =  solution.x - old_solution.x;
        diff_phi      =  solution.phi - old_solution.phi;
        diff_fields   =  [diff_x;diff_phi];
    case {'electro_incompressible','electro_mixed_incompressible'}
        diff_x        =  solution.x.Eulerian_x(:) - old_solution.x.Eulerian_x(:);
        diff_phi      =  solution.phi - old_solution.phi;
        diff_p        =  solution.p - old_solution.p;
        diff_fields   =  [diff_x;diff_phi;diff_p];
    case 'electro_BEM_FEM'
        diff_x        =  solution.x - old_solution.x;
        diff_phi      =  solution.phi - old_solution.phi;
        diff_q0       =  solution.q0 - old_solution.q0;        
        diff_fields   =  [diff_x;diff_phi;diff_q0];
    case {'electro_incompressible_BEM_FEM','electro_mixed_incompressible_BEM_FEM'}
        diff_x        =  solution.x - old_solution.x;
        diff_phi      =  solution.phi - old_solution.phi;
        diff_p        =  solution.p - old_solution.p;
        diff_q0       =  solution.q0 - old_solution.q0;                
        diff_fields   =  [diff_x;diff_phi;diff_p;diff_q0];
end