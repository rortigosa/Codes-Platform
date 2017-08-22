
function material_information     =  GetDerivativesModelMechanics(ielem,dim,ngauss,F,H,J,material_information)                                 

%--------------------------------------------------------------------------
% We get DU: derivatives of the internal energy expressed as a function of
% (F,H,J,D0,d)
%--------------------------------------------------------------------------
material_information                =  FirstDerivativeComputationMechanics(ielem,ngauss,F,H,J,material_information);                                                                      
%--------------------------------------------------------------------------
% We get D2U: second derivatives of the internal energy expressed as a 
% function of (F,H,J,D0,d)
%--------------------------------------------------------------------------
material_information                =  HessianOperatorInternalEnergyMechanics(ielem,dim,ngauss,F,H,J,material_information);
