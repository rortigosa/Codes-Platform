%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function is responsible for the treatment of the contraints in the
% problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str                          =  DirichletBcs(str)

%--------------------------------------------------------------------------
% Constraints for x
%--------------------------------------------------------------------------
switch str.data.formulation
    case 'electro_BEM_FEM'
n_dofs                                =  0;        
    otherwise
[fixdof_u,freedof_u,...
    cons_val_u]                       =  MechanicalDirichletConstraints(str);                         
str.bc.Dirichlet.x.freedof            =  freedof_u;
str.bc.Dirichlet.x.fixdof             =  fixdof_u;
str.bc.Dirichlet.x.cons_val           =  cons_val_u;
n_dofs                                =  size(str.solution.x.Eulerian_x(:),1);
end
%--------------------------------------------------------------------------
% Constraints for phi
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'electro_mechanics','electro_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz','electro_mechanics_Helmholtz_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM',...
          'electro_BEM_FEM'}
[fixdof_phi,freedof_phi,...
    cons_val_phi]                     =  ElectricDirichletConstraints(str);                         
str.bc.Dirichlet.phi.freedof          =  freedof_phi;
str.bc.Dirichlet.phi.fixdof           =  fixdof_phi;
str.bc.Dirichlet.phi.cons_val         =  cons_val_phi;
fixdof_phi                            =  fixdof_phi + n_dofs;
freedof_phi                           =  freedof_phi + n_dofs;
n_dofs                                =  n_dofs + size(str.solution.phi,1);
end
%--------------------------------------------------------------------------
% Constraints for p
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'up','electro_mechanics_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
          freedof_p                   =  (1:size(str.solution.pressure,1))' + n_dofs;
          fixdof_p                    =  [];
          n_dofs                      =  n_dofs + size(str.solution.pressure,1);
end
%--------------------------------------------------------------------------
% Constraints for q0 in BEM/FEM problems
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_BEM_FEM','electro_mechanics_Helmholtz_incompressible_BEM_FEM',...
          'electro_BEM_FEM'}
         freedof_q                    =  n_dofs + (1:size(str.solution.q,1))';                
         fixdof_q                     =  [];
end
%--------------------------------------------------------------------------
% fixdof, freedof and fixed value 
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'u','FHJ','CGC','CGCCascade'}
         str.bc.Dirichlet.fixdof      =  fixdof_u;
         str.bc.Dirichlet.freedof     =  freedof_u;
         str.bc.Dirichlet.cons_val    =  cons_val_u;               
    case 'up'
         str.bc.Dirichlet.fixdof      =  [fixdof_u;   fixdof_p];
         str.bc.Dirichlet.freedof     =  [freedof_u;  freedof_p];
         str.bc.Dirichlet.cons_val    =  cons_val_u;               
    case {'electro_mechanics','electro_mechanics_Helmholtz'}
         str.bc.Dirichlet.fixdof      =  [fixdof_u;    fixdof_phi];
         str.bc.Dirichlet.freedof     =  [freedof_u;   freedof_phi];
         str.bc.Dirichlet.cons_val    =  [cons_val_u;  cons_val_phi];       
    case {'electro_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_Helmholtz_incompressible'}
         str.bc.Dirichlet.fixdof      =  [fixdof_u;    fixdof_phi;   fixdof_p];
         str.bc.Dirichlet.freedof     =  [freedof_u;   freedof_phi;  freedof_p];
         str.bc.Dirichlet.cons_val    =  [cons_val_u;  cons_val_phi];               
    case {'electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz_BEM_FEM'}
         str.bc.Dirichlet.fixdof      =  [fixdof_u;    fixdof_phi;   fixdof_q];
         str.bc.Dirichlet.freedof     =  [freedof_u;   freedof_phi;  freedof_q];
         str.bc.Dirichlet.cons_val    =  [cons_val_u;  cons_val_phi];               
    case {'electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
         str.bc.Dirichlet.fixdof      =  [fixdof_u;    fixdof_phi;   fixdof_p;   fixdof_q];
         str.bc.Dirichlet.freedof     =  [freedof_u;   freedof_phi;  freedof_p;  freedof_q];
         str.bc.Dirichlet.cons_val    =  [cons_val_u;  cons_val_phi];               
    case {'electro_BEM_FEM'}
         str.bc.Dirichlet.freedof     =  [freedof_phi;  freedof_q];
         str.bc.Dirichlet.cons_val    =  cons_val_phi;               
end
end        


