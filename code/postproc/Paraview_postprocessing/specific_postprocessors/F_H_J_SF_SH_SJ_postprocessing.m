function [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff]   =  F_H_J_SF_SH_SJ_postprocessing(gauss_level_information,str)

F                                                 =  gauss_level_information.F;
H                                                 =  gauss_level_information.H;
J                                                 =  gauss_level_information.J';
SigmaF                                            =  gauss_level_information.SigmaF;
SigmaH                                            =  gauss_level_information.SigmaH;
SigmaJ                                            =  gauss_level_information.SigmaJ';
Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
Kirchhoff                                         =  Piola;
