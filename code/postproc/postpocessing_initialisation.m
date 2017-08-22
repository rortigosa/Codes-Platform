function str                                  =  postpocessing_initialisation(str)

dim                                           =  str.data.dim;
nnode                                         =  str.postproc.n_nodes;
%-----------------------------------------
% Eulerian coordinates and electric potential.
%----------------------------------------- 
str.postproc.Eulerian_x                       =  zeros(dim,nnode);
str.postproc.velocity                         =  zeros(dim,nnode);
str.postproc.acceleration                     =  zeros(dim,nnode);
str.postproc.phi                              =  zeros(nnode,1);
%-----------------------------------------
% Second Piola and the different contributions.
%-----------------------------------------
str.postproc.stress                           =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_mec_iso                   =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_mec_aniso                 =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_mec                       =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_diel_iso                  =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_diel_aniso                =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_diel                      =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_piezo                     =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_elec                      =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_deviatoric                =  zeros((dim+1)*dim/2,nnode);
str.postproc.stress_volumetric                =  zeros((dim+1)*dim/2,nnode);
%-----------------------------------------
% Cauchy stress and the different contributions.
%-----------------------------------------
str.postproc.sigma                            =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_mec_iso                    =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_mec_aniso                  =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_mec                        =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_diel_iso                   =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_diel_aniso                 =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_diel                       =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_piezo                      =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_elec                       =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_deviatoric                 =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_volumetric                 =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_pressure                   =  zeros(nnode,1);
%-----------------------------------------
% SigmaF, SigmaH, SigmaJ.
%-----------------------------------------
str.postproc.SigmaF                           =  zeros(dim,dim,nnode);
str.postproc.SigmaH                           =  zeros(dim,dim,nnode);
str.postproc.SigmaJ                           =  zeros(nnode,1);
str.postproc.First_Piola                      =  zeros(dim,dim,nnode);
%-----------------------------------------
% Eulerian and Lagrangian electric field.
%-----------------------------------------
str.postproc.E0                               =  zeros(dim,nnode);
str.postproc.E                                =  zeros(dim,nnode);
str.postproc.E_norm                           =  zeros(nnode,1);
str.postproc.electric_breakdown_factor        =  zeros(nnode,1);
%-----------------------------------------
% Eulerian and Lagrangian electric displacement.
%-----------------------------------------
str.postproc.D                                =  zeros(dim,nnode);
str.postproc.D_diel_iso                       =  zeros(dim,nnode);
str.postproc.D_diel_aniso                     =  zeros(dim,nnode);
str.postproc.D_diel                           =  zeros(dim,nnode);
str.postproc.D_piezo                          =  zeros(dim,nnode);
str.postproc.D0                               =  zeros(dim,nnode);
str.postproc.D0_diel_iso                      =  zeros(dim,nnode);
str.postproc.D0_diel_aniso                    =  zeros(dim,nnode);
str.postproc.D0_diel                          =  zeros(dim,nnode);
str.postproc.D0_piezo                         =  zeros(dim,nnode);
%-----------------------------------------
% Deformation gradient, cofactor and jacobian.
%-----------------------------------------
str.postproc.F                                =  zeros(dim,str.data.dim,nnode);
str.postproc.H                                =  zeros(dim,str.data.dim,nnode);
str.postproc.J                                =  zeros(nnode,1);
str.postproc.FL2norm                          =  zeros(nnode,1);
str.postproc.gradJ                            =  zeros(dim,nnode);
str.postproc.ratioJ                           =  zeros(nnode,1);
%-----------------------------------------
% Energy density function.
%-----------------------------------------
str.postproc.Phi                              =  zeros(nnode,1);
str.postproc.Phi_mec_iso                      =  zeros(nnode,1);
str.postproc.Phi_mec_aniso                    =  zeros(nnode,1);
str.postproc.Phi_mec                          =  zeros(nnode,1);
str.postproc.Phi_diel_iso                     =  zeros(nnode,1);
str.postproc.Phi_diel_aniso                   =  zeros(nnode,1);
str.postproc.Phi_diel                         =  zeros(nnode,1);
str.postproc.Phi_piezo                        =  zeros(nnode,1);
str.postproc.Phi_elec                         =  zeros(nnode,1);
%-----------------------------------------
% Electric density vector.
%-----------------------------------------
str.postproc.J_electric                       =  zeros(dim,nnode);
str.postproc.J0_electric                      =  zeros(dim,nnode);
%-----------------------------------------
% Eigenvalues of the second derivative of
% the internal energy with respect to the
% deformation gradient F.
%-----------------------------------------
str.postproc.eigenvalues                      =  zeros(dim^2,nnode);
str.postproc.E                                =  zeros(dim,nnode);
