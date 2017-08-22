%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function computes:
% The first derivative of the internal energy for electromechanical
% problems in terms of F, H, J, D0 and d
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function  material_information   =  FirstDerivativeComputationMechanics(ielem,n_gauss,F,H,J,material_information)    
%--------------------------------------------------------------------------
% Different models
%--------------------------------------------------------------------------
switch material_information.material_model{material_information.material_identifier(ielem)}
    %----------------------------------------------------------------------
    % Simplified ideal dielectric elastomer
    %----------------------------------------------------------------------
    case 'compressible_Mooney_Rivlin'
        material_information     =  FirstDerivativesCompressibleMooneyRivlin(ielem,n_gauss,...
                                                                  F,H,J,material_information);
end        
        
        