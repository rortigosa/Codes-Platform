%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Specify nodal loads
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str     =  NeumannBcs(str)

switch str.data.formulation
    case {'u','up','FHJ','CGC','CGCCascade'}
         str     =  NeumannBcsMechanics(str);        
    case {'electro_mechanics','electro_mechanics_incompressible','electro_mechanics_mixed_incompressible',...
          'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz','electro_mechanics_Helmholtz_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
         str     =  NeumannBcsElectroMechanics(str);              
    case {'electro_BEM_FEM'}
         str     =  NeumannBcsElectro(str);        
end
