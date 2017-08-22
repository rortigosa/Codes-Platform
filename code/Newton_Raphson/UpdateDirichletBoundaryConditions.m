%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function updates the constrained variables for the next load
% increment.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function  str     =   UpdateDirichletBoundaryConditions(str)

%--------------------------------------------------------------------------
% Update Dirichlet boundary conditions
%--------------------------------------------------------------------------
str.solution.x.Eulerian_x(str.bc.Dirichlet.x.fixdof)   =  str.solution.x.Eulerian_x(str.bc.Dirichlet.x.fixdof) + str.NR.load_factor*str.bc.Dirichlet.x.cons_val;
str.solution.phi(str.bc.Dirichlet.phi.fixdof)          =  str.NR.accumulated_factor*str.bc.Dirichlet.phi.cons_val;
