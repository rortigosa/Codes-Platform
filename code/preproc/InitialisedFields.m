%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Store Lagrangian and Eulerian coordinates. Add contribution from
% Dirichlet BCs just to the Eulerian coordinates. The Lagrangian
% coordinates remain untouched.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                          =  InitialisedFields(str) 
  
switch str.data.formulation
    case 'u'
        str.solution                  =  InitialisedVariablesU(str.geometry,str.mesh);        
    case 'up'
        str.solution                  =  InitialisedVariablesUP(str.geometry,str.mesh);        
    case 'FHJ'
        str.solution                  =  InitialisedVariablesFHJ(str.geometry,str.mesh,str.material_information);        
    case {'CGC','CGCCascade'}
        str.solution                  =  InitialisedVariablesCGC(str.geometry,str.mesh,str.material_information);        
    case {'electro_mechanics','electro_mechanics_Helmholtz'}
        str.solution                  =  InitialisedVariablesElectro(str.geometry,str.mesh,str.quadrature);
    case {'electro_mechanics_incompressible','electro_mechanics_Helmholtz_incompressible'}
        str.solution                  =  InitialisedVariablesElectroIncompressible(str.geometry,str.mesh,str.quadrature);
    case 'electro_mechanics_mixed_incompressible'
        str.solution                  =  InitialisedVariablesElectroMixedIncompressible(str.geometry,str.mesh,str.material_information);
    case {'electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz_BEM_FEM'}
        str.solution                  =  InitialisedVariablesElectroBEMFEM(str.geometry,str.mesh,str.quadrature);
    case {'electro_mechanics_incompressible_BEM_FEM','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
        str.solution                  =  InitialisedVariablesElectroIncompressibleBEMFEM(str.geometry,str.mesh,str.quadrature);
    case 'electro_mechanics_mixed_incompressible_BEM_FEM'
        str.solution                  =  InitialisedVariablesElectroMixedIncompressibleBEMFEM(str.geometry,str.mesh,str.material_information);
    case {'electro_BEM_FEM'}
        str.solution                  =  InitialisedVariablesElectroBEMOnlyElectro(str.geometry,str.mesh,str.quadrature);
    case {'electro'}
        str.solution                  =  InitialisedVariablesOnlyElectro(str.geometry,str.mesh,str.quadrature);
end                                  
%--------------------------------------------------------------------------
% Initialisation of the incremental solution for the Newton-Raphson algorithm
%--------------------------------------------------------------------------
str.solution.incremental_solution     =  zeros(str.solution.n_dofs,1);
%--------------------------------------------------------------------------
%  Lagrange multipliers for contact
%--------------------------------------------------------------------------
if str.contact.lagrange_multiplier
   str.solution.contact_multiplier    =  zeros(str.mesh.surface.contact_multiplier.n_nodes,1);
end



