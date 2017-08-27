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
         str.solution         =  UpdatedIncrementalVariablesU(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case 'FHJ'
         str.solution         =  UpdatedIncrementalVariablesFHJ(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case 'CGC'
         str.solution         =  UpdatedIncrementalVariablesCGC(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case 'CGCCascade'
         str.solution         =  UpdatedIncrementalVariablesCGCCascade(str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro'}
         str.solution         =  UpdatedIncrementalVariablesOnlyElectro(str.solution);
    case {'electro_BEM_FEM'}
         str.solution         =  UpdatedIncrementalVariablesElectroBEMFEMOnlyElectro(str.mesh,str.solution);
    case {'electro_mechanics','electro_mechanics_Helmholtz'}
         str.solution         =  UpdatedIncrementalVariablesElectro(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro_mechanics_incompressible','electro_mechanics_Helmholtz_incompressible'}
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressible(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mechanics_mixed_incompressible'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompressible(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case {'electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz_BEM_FEM'}
         str.solution         =  UpdatedIncrementalVariablesElectroBEMFEM(str.geometry.dim,str.mesh,str.solution);
    case 'electro_mechanics_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressibleBEMFEM(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mechanics_mixed_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompBEMFEM(str.geometry.dim,str.mesh,str.solution,str.time_integrator,str.contact);
end


