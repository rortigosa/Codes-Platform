function [I,W_dot_elec]   =  analytical_Electric_Intensity(str)


%--------------------------------------------------------------------------
% Material properties.
%--------------------------------------------------------------------------
vacuum_permittivity       =  8.854187817620e-12;
c1                        =  max(max(str.properties.im1));
c2                        =  max(max(str.properties.im2));
c3                        =  max(max(str.properties.im3));
c12                       =  max(c1,c2);
c13                       =  max(c1,c3);
cf                        =  max(c12,c13);
c                         =  1/(cf*sqrt(vacuum_permittivity));
lambda_iso                =  str.properties.lambda_iso(1,1)*c^2;
%--------------------------------------------------------------------------
% Equivalent electric permittivity.
%--------------------------------------------------------------------------
epsilon                   =  2*lambda_iso;
%--------------------------------------------------------------------------
% Area of the capacitor.
%--------------------------------------------------------------------------
Xmax                      =  max(str.Lagrangian_X(1,:));
Xmin                      =  min(str.Lagrangian_X(1,:));
A                         =  Xmax - Xmin;      
%--------------------------------------------------------------------------
% Lenth of the capacitor.
%--------------------------------------------------------------------------
Ymax                      =  max(str.Lagrangian_X(2,:));
Ymin                      =  min(str.Lagrangian_X(2,:));
L                         =  Ymax - Ymin;
%--------------------------------------------------------------------------
% Capacitance of the capacitor.
%--------------------------------------------------------------------------
C                         =  epsilon*A/L;
%--------------------------------------------------------------------------
% DphiDt.
%--------------------------------------------------------------------------
phi                       =  str.phi(str.plotting.phi_node,1);
phi_n1                    =  str.phin1(str.plotting.phi_node,1);
Dt                        =  str.data.delta_t;
DphiDt                    =  (phi - phi_n1)/(Dt);
%--------------------------------------------------------------------------
% Intensity and electric power.
%--------------------------------------------------------------------------
I                         =  C*DphiDt;
W_dot_elec                =  -I*phi;

