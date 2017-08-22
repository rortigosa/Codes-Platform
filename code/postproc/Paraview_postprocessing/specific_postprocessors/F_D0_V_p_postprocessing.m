function [option,F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff,D0,E0,pressure]  =  F_D0_V_p_postprocessing(str,gauss_level_information,E0)

F                                                 =  gauss_level_information.F;
H                                                 =  gauss_level_information.H;
J                                                 =  gauss_level_information.J';
SigmaF                                            =  gauss_level_information.SigmaF;
SigmaH                                            =  gauss_level_information.SigmaH;
SigmaJ                                            =  gauss_level_information.SigmaJ';
Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
Kirchhoff                                         =  Piola;
SigmaD0                                           =  gauss_level_information.DUDD0;
SigmaV                                            =  gauss_level_information.DUDV;
D0                                                =  gauss_level_information.D0;
V                                                 =  gauss_level_information.V;
pressure                                          =  (gauss_level_information.P)';
for igauss=1:size(str.quadrature.Chi,1)
    E0(:,igauss)                                  =  SigmaD0(:,igauss) + F(:,:,igauss)'*SigmaV(:,igauss);
end
option  =  cell(1,3);
option{1,1}                                       =  SigmaV;
option{1,2}                                       =  D0;
option{1,3}                                       =  V;
