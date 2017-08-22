function [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff,E0,D0] =  F_D0_postprocessing(str,gauss_level_information)

F                                                 =  gauss_level_information.F;
H                                                 =  gauss_level_information.H;
J                                                 =  gauss_level_information.J';
SigmaF                                            =  gauss_level_information.SigmaF;
SigmaH                                            =  gauss_level_information.SigmaH;
SigmaJ                                            =  gauss_level_information.SigmaJ';
Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
Kirchhoff                                         =  Piola;
E0                                                =  gauss_level_information.DUDD0;
D0                                                =  gauss_level_information.D0;
