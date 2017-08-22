%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cauchy stress for postprocessing. We have to be very carefull because the
% resulting Cauchy stress needs to satisfy material frame indefference.
% Therefore, it needs to be symmetric.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma]                                  = cauchy_stress(str,Kirchhoff,F_,H_,J_,SigmaF_,SigmaH_,SigmaJ_,optional,gauss_level_information)

switch str.data.formulation_type
        case {'displacement_potential_formulation','u_p_phi_compressible'}
             sigma                                =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
             for igauss=1:size(str.postproc.f_e.N,2)
                 F                                =  str.grad.F(:,:,igauss);
                 J                                =  det(F);
                 sigma(:,:,igauss)                =  1/J*F*Kirchhoff(:,:,igauss)*F';
             end
        case 'mixed_formulation'
             sigma                                =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
%            sigma_                               =  sigma;
             switch str.data.formulation.mixed_type
                 case {'full_mixed_formulation','H_J_mixed_formulation','complementary_energy','stabilised_complementary_energy',...
                       'stabilised_complementary_energy_1','stabilised_complementary_energy_2','stabilised_complementary_energy_3','displacement_formulation','displacement_potential_mixed_formulation',...
                       'full_mixed_formulation_electroelasticity_F_D0','full_mixed_formulation_electroelasticity_F_E0',...
                       'full_mixed_formulation_electroelasticity_S_D0','full_mixed_formulation_electroelasticity_S_E0','u_p_phi_incompressible'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          H                       =  H_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaH                  =  SigmaH_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauH                    =  Javier_double_cross_product(SigmaH*H',eye(str.data.dim),1,1,str.data.dim);
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)       =  1/J*(tauF + tauH + tauJ);
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
%                         sigma(:,:,igauss)       =  1/J*Kirchhoff(:,:,igauss)*F';
                      end                      
                 case {'F_J_mixed_formulation'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)       =  1/J*(tauF + tauJ);
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
%                         sigma(:,:,igauss)       =  1/J*Kirchhoff(:,:,igauss)*F';
                      end                      
                 case {'u_J_mixed_formulation'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)       =  1/J*(tauF + tauJ);
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
%                         sigma(:,:,igauss)       =  1/J*Kirchhoff(:,:,igauss)*F';
                      end                      
                 case {'full_mixed_formulation_electroelasticity_F_D0_V','full_mixed_formulation_electroelasticity_F_E0_V'}
                          SigmaV_                 =  optional{1,1};
                          D0_                     =  optional{1,2};
                          V_                      =  optional{1,3};
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          H                       =  H_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaV                  =  SigmaV_(:,igauss);
                          D0                      =  D0_(:,igauss);
                          V                       =  V_(:,igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaH                  =  SigmaH_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauH                    =  Javier_double_cross_product(SigmaH*H',eye(str.data.dim),1,1,str.data.dim);
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)        =  1/J*(tauF + tauH + tauJ + SigmaV*V');
                         
                         
%                           SigmaF                  =  SigmaF_(:,:,igauss);
%                           SigmaH                  =  SigmaH_(:,:,igauss);
%                           SigmaJ                  =  SigmaJ_(igauss);                          
%                           tauF                    =  SigmaF*str.grad.F(:,:,igauss)';
%                           tauH                    =  Javier_double_cross_product(SigmaH*str.grad.H(:,:,igauss)',eye(str.data.dim),1,1,str.data.dim);
%                           tauJ                    =  str.grad.J(igauss)*SigmaJ*eye(str.data.dim);
%                          sigma(:,:,igauss)        =  1/str.grad.J(igauss)*(tauF + tauH + tauJ + SigmaV*(str.grad.F(:,:,igauss)*D0)');
                         
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
%                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
                         %sigma(:,:,igauss)       =  1/J*Kirchhoff(:,:,igauss)*F';
                      end                      
                 case {'full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          H                       =  H_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaH                  =  SigmaH_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauH                    =  Javier_double_cross_product(SigmaH*H',eye(str.data.dim),1,1,str.data.dim);
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)       =  1/J*(tauF + tauH + tauJ) + gauss_level_information.p(igauss)*eye(str.data.dim);
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
%                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
                      end                      
                 case {'full_incompressible','nearly_incompressible_three_field','nearly_incompressible_two_field','compressible_two_field'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F3D                     =  F_(:,:,igauss);
                      switch str.data.dim
                             case 2
                                  F3D(3,3)        =  1;
                      end
                      J                           =  J_(igauss);
                      sigma(:,:,igauss)           =  1/J*Kirchhoff(:,:,igauss)*F3D';
                      end
                 case {'twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
                      for igauss=1:size(str.postproc.f_e.N,2)
                          F                       =  F_(:,:,igauss);
                          H                       =  H_(:,:,igauss);
                          J                       =  J_(igauss);
                          SigmaF                  =  SigmaF_(:,:,igauss);
                          SigmaJ                  =  SigmaJ_(igauss);                          
                          tauF                    =  SigmaF*F';
                          tauJ                    =  J*SigmaJ*eye(str.data.dim);
                         sigma(:,:,igauss)       =  1/J*(tauF + tauJ);
                          FF                      =  str.grad.F(:,:,igauss);
                          JJ                      =  str.grad.J(igauss);
                         sigma(:,:,igauss)       =  1/JJ*Kirchhoff(:,:,igauss)*FF';
                      end
             end
end
