%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% The fields of a specific formulation are updated here
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                  =  FieldsUpdate(str)

switch str.data.formulation
    case 'u'
         str.solution         =  UpdatedIncrementalVariablesU(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'FHJ'
         str.solution         =  UpdatedIncrementalVariablesFHJ(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'CGC'
         str.solution         =  UpdatedIncrementalVariablesCGC(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'CGCCascade'
         str.solution         =  UpdatedIncrementalVariablesCGCCascade(str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro'}
         str.solution         =  UpdatedIncrementalVariablesOnlyElectro(str.solution);
    case {'electro_BEM_FEM'}
         str.solution         =  UpdatedIncrementalVariablesElectroBEMFEMOnlyElectro(str.mesh,str.solution);
    case {'electro_mechanics','electro_mechanics_Helmhotlz'}
         str.solution         =  UpdatedIncrementalVariablesElectro(str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro_mechanics_incompressible','electro_mechanics_Helmholtz_incompressible'}
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressible(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mechanics_mixed_incompressible'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompressible(str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz_BEM_FEM'}
         str.solution         =  UpdatedIncrementalVariablesElectroBEMFEM(str.mesh,str.solution);
    case 'electro_mechanics_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressibleBEMFEM(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mechanics_mixed_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompBEMFEM(str.mesh,str.solution,str.time_integrator,str.contact);
end


