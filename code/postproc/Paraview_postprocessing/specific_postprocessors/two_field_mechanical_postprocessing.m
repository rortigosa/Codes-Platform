function [F,H,J,Kirchhoff,p]                      =  two_field_mechanical_postprocessing(str)

str.f_e.N                                         =  str.postproc.f_e.N;
str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
str                                               =  gradients(xelem,Xelem,phielem,str);
F                                                 =  str.grad.F;
H                                                 =  str.grad.H;
J                                                 =  str.grad.J;
PPiola                                            =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
for igauss=1:size(str.quadrature.Chi,1)
    Piola(:,:,igauss)                             =  PPiola(:,:,igauss) + 0*gauss_level_information.p(igauss)*eye(3);
end
Kirchhoff                                         =  Piola;
p                                                 =  gauss_level_information.p;
