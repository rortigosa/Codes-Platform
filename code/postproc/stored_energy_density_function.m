%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the stored energy function for the different
% material models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function   [str]      =  stored_energy_density_function(str)

%--------------------------------------------------------------------------
% Selecting the material model
%--------------------------------------------------------------------------
filename              =  'material_model.rst'; 
fid                   =  fopen(filename, 'r');  
material_model        =  fgetl(fid);
fclose(fid);
 
switch material_model
    case 'schroder'
         [str]        =  stored_energy_density_schroder(str);          
    case 'paper_form_transv_iso_rhombic'
         [str]        =  stored_energy_density_transverse_paper_form(str);
    case 'paper_form_isotropy'
         [str]        =  stored_energy_density_isotropic_paper_form(str);          
end

