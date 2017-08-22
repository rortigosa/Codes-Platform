function [F,H,J,E0,D0,Kirchhoff]        =  displacement_potential_postprocessing(str,ielem)

oldstr                                                        =  str;
str.f_e.N                                                     =  str.postproc.f_e.N;
str.f_e.DN_chi                                                =  str.postproc.f_e.DN_chi;
new_str                                                       =  gradients(xelem,Xelem,phielem,str);
str                                                           =  oldstr;
J                                                             =  new_str.grad.J;
F                                                             =  new_str.grad.F;
H                                                             =  new_str.grad.H;
E0                                                            =  new_str.grad.E0;
E                                                             =  new_str.grad.E;
%-----------------------------------------
% Second Piola and different component.
%-----------------------------------------
switch str.vacuum.flag
    case 0
        switch str.formulation
            case 'D_to_E'
                [S]                                           =  second_piola_d(str);
            case 'E'
                [S]                                           =  second_piola(str);
        end
    case 1
        vacuum_identifier                                     =  max(str.properties.material);
        switch str.properties.material(ielem)
            case vacuum_identifier
                str.formulation                               =  'E';
                [S]                                           =  second_piola_vacuum(str);
            otherwise
                str.formulation                               =  str.old_formulation;
                [S]                                           =  second_piola_d(str);
        end
end
Kirchhoff                                                     =  S;
%-----------------------------------------
% Lagrangian electric displacement.
%-----------------------------------------
str.ielem                                                     =  ielem;
switch str.vacuum.flag
    case 0
        switch str.formulation
            case 'D_to_E'
                [D0,str]                                      =  internal_Newton_Raphson(str);
            case 'E'
                D0                                            =  electric_displacement(str);
        end
    case 1
        vacuum_identifier                                     =  max(str.properties.material);
        switch str.properties.material(ielem)
            case vacuum_identifier
                str.formulation                               =  'E';
                [D0]                                          =  electric_displacement_vacuum(str);
            otherwise
                str.formulation                               =  str.old_formulation;
                [D0,str]                                      =  internal_Newton_Raphson(str);
        end
end
