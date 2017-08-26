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
switch str.data.formulation
    case {'u','up','FHJ','CGC','CGCCascade','electro_mechanics','electro_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz','electro_mechanics_Helmholtz_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
str.solution.x.Eulerian_x(str.bc.Dirichlet.x.fixdof)   =  str.solution.x.Eulerian_x(str.bc.Dirichlet.x.fixdof) + str.NR.load_factor*str.bc.Dirichlet.x.cons_val;
end
switch str.data.formulation
    case {'electro_mechanics','electro_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz','electro_mechanics_Helmholtz_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM',...
          'electro_BEM_FEM','electro'}
str.solution.phi(str.bc.Dirichlet.phi.fixdof)          =  str.NR.accumulated_factor*str.bc.Dirichlet.phi.cons_val;
end      
      
      