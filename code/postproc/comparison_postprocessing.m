%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is needed when two examples are going to be compared. This 
% is the postprocessing of variables needed to plot. This function is a
% postprocessing for the error or difference between the two methods/models
% for the same example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [str]                                      =  comparison_postprocessing(str)
   
str.postproc.displacement                           =  variable_error(str.master.Eulerian_x,str.slave.Eulerian_x);
str.postproc.phi                                    =  variable_error(str.master.phi,str.slave.phi);
str.postproc.stress                                 =  variable_error(str.master_post.stress,str.slave_post.stress);
str.postproc.stress_mec_iso                         =  variable_error(str.master_post.stress_mec_iso,str.slave_post.stress_mec_iso);
str.postproc.stress_mec_aniso                       =  variable_error(str.master_post.stress_mec_aniso,str.slave_post.stress_mec_aniso);
str.postproc.stress_mec                             =  variable_error(str.master_post.stress_mec,str.slave_post.stress_mec);
str.postproc.stress_diel_iso                        =  variable_error(str.master_post.stress_diel_iso,str.slave_post.stress_diel_iso);
str.postproc.stress_diel_aniso                      =  variable_error(str.master_post.stress_diel_aniso,str.slave_post.stress_diel_aniso);
str.postproc.stress_diel                            =  variable_error(str.master_post.stress_diel,str.slave_post.stress_diel);
str.postproc.stress_piezo                           =  variable_error(str.master_post.stress_piezo,str.slave_post.stress_piezo);
str.postproc.stress_elec                            =  variable_error(str.master_post.stress_elec,str.slave_post.stress_elec);
str.postproc.stress_deviatoric                      =  variable_error(str.master_post.stress_deviatoric,str.slave_post.stress_deviatoric);
str.postproc.stress_volumetric                      =  variable_error(str.master_post.stress_volumetric,str.slave_post.stress_volumetric);
str.postproc.sigma                                  =  variable_error(str.master_post.sigma,str.slave_post.sigma);
str.postproc.sigma_mec_iso                          =  variable_error(str.master_post.sigma_mec_iso,str.slave_post.sigma_mec_iso);
str.postproc.sigma_mec_aniso                        =  variable_error(str.master_post.sigma_mec_aniso,str.slave_post.sigma_mec_aniso);
str.postproc.sigma_mec                              =  variable_error(str.master_post.sigma_mec,str.slave_post.sigma_mec);
str.postproc.sigma_diel_iso                         =  variable_error(str.master_post.sigma_diel_iso,str.slave_post.sigma_diel_iso);
str.postproc.sigma_diel_aniso                       =  variable_error(str.master_post.sigma_diel_aniso,str.slave_post.sigma_diel_aniso);
str.postproc.sigma_diel                             =  variable_error(str.master_post.sigma_diel,str.slave_post.sigma_diel);
str.postproc.sigma_piezo                            =  variable_error(str.master_post.sigma_piezo,str.slave_post.sigma_piezo);
str.postproc.sigma_elec                             =  variable_error(str.master_post.sigma_elec,str.slave_post.sigma_elec);
str.postproc.sigma_deviatoric                       =  variable_error(str.master_post.sigma_deviatoric,str.slave_post.sigma_deviatoric);
str.postproc.sigma_volumetric                       =  variable_error(str.master_post.sigma_volumetric,str.slave_post.sigma_volumetric);
str.postproc.sigma_pressure                         =  variable_error(str.master_post.sigma_pressure,str.slave_post.sigma_pressure);
str.postproc.E0                                     =  variable_error(str.master_post.E0,str.slave_post.E0);
str.postproc.E                                      =  variable_error(str.master_post.E,str.slave_post.E);
str.postproc.D                                      =  variable_error(str.master_post.D,str.slave_post.D);
str.postproc.D_diel_iso                             =  variable_error(str.master_post.D_diel_iso,str.slave_post.D_diel_iso);
str.postproc.D_diel_aniso                           =  variable_error(str.master_post.D_diel_aniso,str.slave_post.D_diel_aniso);
str.postproc.D_diel                                 =  variable_error(str.master_post.D_diel,str.slave_post.D_diel);
str.postproc.D_piezo                                =  variable_error(str.master_post.D_piezo,str.slave_post.D_piezo);
str.postproc.D0                                     =  variable_error(str.master_post.D0,str.slave_post.D0);
str.postproc.D0_diel_iso                            =  variable_error(str.master_post.D0_diel_iso,str.slave_post.D0_diel_iso);
str.postproc.D0_diel_aniso                          =  variable_error(str.master_post.D0_diel_aniso,str.slave_post.D0_diel_aniso);
str.postproc.D0_diel                                =  variable_error(str.master_post.D0_diel,str.slave_post.D0_diel);
str.postproc.D0_piezo                               =  variable_error(str.master_post.D0_piezo,str.slave_post.D0_piezo);
str.postproc.F                                      =  variable_error(str.master_post.F,str.slave_post.F);
str.postproc.J                                      =  variable_error(str.master_post.J,str.slave_post.J);
str.postproc.Phi                                    =  variable_error(str.master_post.Phi,str.slave_post.Phi);
str.postproc.Phi_mec_iso                            =  variable_error(str.master_post.Phi_mec_iso,str.slave_post.Phi_mec_iso);
str.postproc.Phi_mec_aniso                          =  variable_error(str.master_post.Phi_mec_aniso,str.slave_post.Phi_mec_aniso);
str.postproc.Phi_mec                                =  variable_error(str.master_post.Phi_mec,str.slave_post.Phi_mec);
str.postproc.Phi_diel_iso                           =  variable_error(str.master_post.Phi_diel_iso,str.slave_post.Phi_diel_iso);
str.postproc.Phi_diel_aniso                         =  variable_error(str.master_post.Phi_diel_aniso,str.slave_post.Phi_diel_aniso);
str.postproc.Phi_diel                               =  variable_error(str.master_post.Phi_diel,str.slave_post.Phi_diel);
str.postproc.Phi_piezo                              =  variable_error(str.master_post.Phi_piezo,str.slave_post.Phi_piezo);
str.postproc.Phi_elec                               =  variable_error(str.master_post.Phi_elec,str.slave_post.Phi_elec);
str.postproc.J_electric                             =  variable_error(str.master_post.J_electric,str.slave_post.J_electric);
str.postproc.J0_electric                            =  variable_error(str.master_post.J0_electric,str.slave_post.J0_electric);
end


function [error]                                    =  variable_error(master_variable,slave_variable)
%   error                                           =  abs(master_variable - slave_variable)./master_variable*100;     
%   [error]                                         =  NaN_Inf_tensor_components(error);    
    denominator                                     =  max(max(max(master_variable)));
    if denominator==0
       denominator                                  =  1;
    end
    error                                           =  abs(master_variable - slave_variable)/denominator*100;     
end









