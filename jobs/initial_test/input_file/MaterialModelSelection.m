%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function selects the constitutive model used for the material.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str                                  =  MaterialModelSelection(str)

str.material_information.n_material_models    =  1;
%material_model{1}                            =  'isochoric_Mooney_Rivlin';
%material_model{2}                            =  'isochoric_Mooney_Rivlin';
%material_model{1}                            =  'compressible_Mooney_Rivlin';
%material_model{1}                            =  'polyconvex_transversely_isotropic';
%material_model{1}                            =  'ideal_dielectric_elastomer';
material_model{1}                             =  'simplified_ideal_dielectric_elastomer';
material_model{1}                             =  'ideal_dielectric_elastomer';

if size(material_model,2)==str.material_information.n_material_models
else
   fprintf('\n Incorrect definition of the material model. Go to function "material_model_selection.m" \n')   
end

str.material_information.material_model       =  material_model;

str.material_information.material_identifier  =  ones(str.mesh.volume.n_elem,1);

if max(str.material_information.material_identifier)==str.material_information.n_material_models
else
   fprintf('\n Incorrect definition of the material model. Go to function "material_model_selection.m" \n')
end
