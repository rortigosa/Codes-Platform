%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Initialisation of the formulation.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str   =  InitialisationFormulation(str)

%--------------------------------------------------------------------------
% Initialisation of the fields associated with the chosen formulation
%--------------------------------------------------------------------------  
str            =  InitialisedFields(str); 
%--------------------------------------------------------------------------
% Initialise constant contributions the mass matrix
%--------------------------------------------------------------------------  
switch str.data.analysis
    case 'dynamic'
switch str.data.formulation
    case {'electro_BEM_FEM','electro'}
    otherwise
str.assembly   =  ConstantMassMatrices(str.geometry,str.mesh,str.fem,str.quadrature);
str            =  MassMatricesAssembly(str);
end
end
%--------------------------------------------------------------------------
% Preprocessing the vectorisation of the code                                                                                                                                                                                                                                                                                                                
%--------------------------------------------------------------------------
str            =  VectorisationMatlabCode(str);
%--------------------------------------------------------------------------
% Initialisation of parameters for the Newton-Raphson algorithm
%--------------------------------------------------------------------------
str            =  NRInitialisation(str);
%--------------------------------------------------------------------------
% Initialisation for contact problems
%--------------------------------------------------------------------------
str            =  ContactInitialisation(str);
%--------------------------------------------------------------------------
% Initialisation for residuals
%--------------------------------------------------------------------------
str            =  InitialisedResiduals(str);



