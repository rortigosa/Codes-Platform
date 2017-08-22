%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% The fields of a specific formulation are updated here
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                  =  FieldsUpdate(str)

switch str.data.formulation
    case 'electro'
         str.solution         =  UpdatedIncrementalVariablesElectro(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_incompressible'
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressible(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mixed_incompressible'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompressible(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroBEMFEM(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroIncompressibleBEMFEM(str.mesh,str.solution,str.time_integrator,str.contact);
    case 'electro_mixed_incompressible_BEM_FEM'
         str.solution         =  UpdatedIncrementalVariablesElectroMixedIncompBEMFEM(str.mesh,str.solution,str.time_integrator,str.contact);
end