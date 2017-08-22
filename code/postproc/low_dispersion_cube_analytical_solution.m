function [U,p,J,Piola]                               =  low_dispersion_cube_analytical_solution(XX,t,str)
% %--------------------------------------------------------------------------
% % Constraint all dofs in X=0
% %--------------------------------------------------------------------------
U                                            =  zeros(3,1);

mu                                           =  6.538461538461538e+06;
lambda                                       =  9.807692307692308e+06;
rho0                                         =  1100;
A                                            =  1;
B                                            =  1;
C                                            =  1;
U0                                           =  5e-4;

X                                        =  XX(1);
Y                                        =  XX(2);
Z                                        =  XX(3);

UX                                       =  A*sin(pi*X/2)*cos(pi*Y/2)*cos(pi*Z/2);
UY                                       =  B*cos(pi*X/2)*sin(pi*Y/2)*cos(pi*Z/2);
UZ                                       =  C*cos(pi*X/2)*cos(pi*Y/2)*sin(pi*Z/2);
 
U                                        =  [UX;UY;UZ];
 
cd                                       =  sqrt((2*mu + lambda)/rho0);
%cd                                       =  sqrt((mu)/rho0);
U                                        =  U0*cos(sqrt(3)/2*cd*pi*t)*U;


F                                        =  zeros(3,3);
F(1,1)                                   =  U0*cos(sqrt(3)/2*cd*pi*t)*A*pi/2*cos(pi*X/2)*cos(pi*Y/2)*cos(pi*Z/2) + 1;
F(1,2)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*A*pi/2*sin(pi*X/2)*sin(pi*Y/2)*cos(pi*Z/2);
F(1,3)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*A*pi/2*sin(pi*X/2)*cos(pi*Y/2)*sin(pi*Z/2);

F(2,1)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*B*pi/2*sin(pi*X/2)*sin(pi*Y/2)*cos(pi*Z/2);
F(2,2)                                   =  U0*cos(sqrt(3)/2*cd*pi*t)*B*pi/2*cos(pi*X/2)*cos(pi*Y/2)*cos(pi*Z/2)+ 1;
F(2,3)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*B*pi/2*cos(pi*X/2)*sin(pi*Y/2)*sin(pi*Z/2);


F(3,1)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*C*pi/2*sin(pi*X/2)*cos(pi*Y/2)*sin(pi*Z/2) ;
F(3,2)                                   =  -U0*cos(sqrt(3)/2*cd*pi*t)*C*pi/2*cos(pi*X/2)*sin(pi*Y/2)*sin(pi*Z/2);
F(3,3)                                   =  U0*cos(sqrt(3)/2*cd*pi*t)*C*pi/2*cos(pi*X/2)*cos(pi*Y/2)*cos(pi*Z/2) + 1;

J                                        =  det(F);


%k                                        =  1.416666666666667e+07;
k                                        =  str.properties.material_parameters.k;
p                                        =  k*(J - 1);


H   =  J*inv(F)';
        eta =  str.properties.material_parameters.eta;
        gamma =  str.properties.material_parameters.gamma;
%        eta =  str.properties.material_parameters.alpha;
%        gamma =  str.properties.material_parameters.beta;
SigmaF                                      =  2*eta/(J^(2/3))*F;
SigmaH                                      =  3*gamma/J^2*(H'*H)^(1/2)*H;
SigmaJ                                      =  -2/3*eta/J^(5/3)*(F'*F) - 2*gamma/J^3*(H'*H)^(3/2);
Piola                                       =  SigmaF + Javier_double_cross_product(SigmaH,F,0,0,3) + (SigmaJ + p)*H;
