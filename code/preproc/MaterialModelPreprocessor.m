%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  This function redas the material model selected and its associated
%  material properties
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function str  =  MaterialModelPreprocessor(str)

%--------------------------------------------------------------------------
% Universal constants
%--------------------------------------------------------------------------
str           =  UniversalPhysicalConstants(str);
%--------------------------------------------------------------------------
% Select the constitutive model
%--------------------------------------------------------------------------
str           =  MaterialModelSelection(str);
%--------------------------------------------------------------------------
% Select the material parameters of the constitutive model
%--------------------------------------------------------------------------
str           =  MaterialParametersInput(str);
%--------------------------------------------------------------------------
% Initialise first and second derivatives of the constitutive model
%--------------------------------------------------------------------------
str           =  InitialisationDerivativesModel(str); 


