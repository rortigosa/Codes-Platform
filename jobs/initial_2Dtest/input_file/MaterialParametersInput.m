
function    str                                                     =  MaterialParametersInput(str)


str.material_information.material_parameters.density(1)             =  1e3;

e0                                                                  =  8.854e-12;
for imat = 1:str.material_information.n_material_models
   switch str.material_information.material_model{imat}    
    case 'compressible_Mooney_Rivlin'
         mu1                                                        =  1e5;
         mu2                                                        =  1e5;
         lambda                                                     =  1e6;
         str.material_information.material_parameters.mu1(imat)     =  mu1;
         str.material_information.material_parameters.mu2(imat)     =  mu2;
         str.material_information.material_parameters.kappa(imat)   =  lambda;
    case 'isochoric_Mooney_Rivlin'
         mu1                                                        =  1e5;
         mu2                                                        =  1e5;
         lambda                                                     =  1e6;
         str.material_information.material_parameters.mu1(imat)     =  mu1;
         str.material_information.material_parameters.mu2(imat)     =  mu2;
         str.material_information.material_parameters.kappa(imat)   =  lambda;
       case 'ideal_dielectric_elastomer'
         mu1                                                        =  1e5;
         mu2                                                        =  1e5;
         lambda                                                     =  1e6;
         e1                                                         =  4*e0;
         e2                                                         =  1e3*e0;
         str.material_information.material_parameters.mu1(imat)     =  mu1;
         str.material_information.material_parameters.mu2(imat)     =  mu2;
         str.material_information.material_parameters.kappa(imat)   =  lambda;
         str.material_information.material_parameters.e1(imat)      =  e1;
         str.material_information.material_parameters.e2(imat)      =  e2;
       case 'simplified_ideal_dielectric_elastomer'
         mu1                                                        =  1e5;
         mu2                                                        =  1e5;
         lambda                                                     =  1e6;
         e1                                                         =  1*e0;
         e2                                                         =  1e3*e0;
         str.material_information.material_parameters.mu1(imat)     =  mu1;
         str.material_information.material_parameters.mu2(imat)     =  mu2;
         str.material_information.material_parameters.kappa(imat)   =  lambda;
         str.material_information.material_parameters.e1(imat)      =  e1;
         str.material_information.material_parameters.e2(imat)      =  e2;           
   end
end

