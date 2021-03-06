%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The influence of all the nodes is considered
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [solution,str]                         =  output_arclength_v5(physics,solution,constraint_dof_u,Rw0_fixed,w0_fixed_residual,p0_fixed_residual,sk,total_u_fixed,str,iteration)


dirichlet  =  1;

arc_length_radius                               =  str.arc_length_radius;

%mechanical_activation                          =  1;
%electric_activation                            =  1;
mechanical_scaling                              =  1;
%potential_scaling                              =  1e5*str.properties.material_parameters.e/sqrt(str.properties.material_parameters.mu*str.properties.material_parameters.e);
potential_scaling                               =  0*1e0*str.properties.material_parameters.e/sqrt(str.properties.material_parameters.mu*str.properties.material_parameters.e);
%scaling                                        =  1/sqrt(str.properties.material_parameters.mu*str.properties.material_parameters.e)/str.data.scaling_area;

scaling                                         =  str.arc_length_scaling;

theta                                           =  str.arc_length.theta;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Stiffness matrices  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
fixed_potential                                 =  zeros(str.n_nodes,1);
fixed_potential(str.elec.freedof_prescribed)    =  str.elec.freedof_prescribed_value;                  
Kxx                                             =  physics.K_matrix;
KxT0                                            =  zeros(size(str.assemb_matrix.total_matrix,1),size(str.freedof_prescribed,1));
KxT0(str.freedof_prescribed,...
    (1:size(str.freedof_prescribed,1))')        =  eye(size(str.freedof_prescribed,1));

KxT0(str.fixdof_arclength,:)                    =  [];
KxT0                                            =  sparse(KxT0);



%--------------------------------------------------------------------------
% Kulambda according to linear Dirichlet boundary conditions
%--------------------------------------------------------------------------
if dirichlet
X_prescribed0                                   =  [str.fixed_displacement(:);fixed_potential];
X_prescribed0_auxiliar                          =  X_prescribed0(str.freedof_prescribed);
X_prescribed0                                   =  X_prescribed0_auxiliar;
Kulambda                                         =  sparse(zeros(size(Kxx,1)+size(KxT0,2),1));
Kulambda(size(Kxx,1)+1:size(Kxx,1)+size(KxT0,2)) =  -X_prescribed0;
else
%--------------------------------------------------------------------------
% Kulambda according to rotational Dirichlet boundary conditions
%--------------------------------------------------------------------------
X_prescribed0                                   =  str.previous_Eulerian_x(str.freedof_prescribed);
X_prescribed0_                                  =  reshape(X_prescribed0,2,[]);
dR                                              =  [[-sin(theta)  -cos(theta)];[cos(theta) -sin(theta)]];
KUlambda                                        =  -dR*X_prescribed0_;
KUlambda                                        =  reshape(KUlambda,[],1);
Kulambda                                        =  sparse(zeros(size(Kxx,1)+size(KxT0,2),1));
Kulambda(size(Kxx,1)+1:size(Kxx,1)+size(KxT0,2)) = KUlambda;
end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Global assembly of the matrix
%--------------------------------------------------------------------------
Kuu                                             =  [Kxx;KxT0'];
rows                                            =  1:size(Kxx,1);
columns                                         =  size(Kxx,1)+1:size(Kxx,1) + size(X_prescribed0,1);
Kuu(rows,columns)                               =  KxT0;

%Ixx                                             =  sparse(zeros(size(Kxx,1),size(Kxx,1)+size(KxT0,2)));
%Ixx(1:size(Kxx,1)+1:size(Kxx,1)^2)              =  ones(size(Kxx,1),1); 

%IT0T0                                           =  sparse(zeros(size(KxT0,2),size(Kxx,1)+size(KxT0,2)));
%initial                                         =  size(KxT0,2)*size(Kxx,1) + 1;
%final                                           =  size(KxT0,2)*(size(Kxx,1) + size(KxT0,2));
%Delta                                           =  size(KxT0,2) + 1;
%IT0T0(initial:Delta:final)                      =  ones(size(KxT0,2),1); 


%invKxx                                          =  Kxx\Identity;
%M                                               =  KxT0'*(Kxx\KxT0);
%Identity_M                                      =  sparse(eye(size(M,1)));
%M                                               =  M\Identity;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Residuals
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Rx                                              =  physics.Residual;

%--------------------------------------------------------------------------
% RT0 according to linear Dirichlet boundary conditions
%--------------------------------------------------------------------------
if dirichlet
RT0                                             =  [str.Eulerian_x(:) - str.Lagrangian_X(:) - str.data.acumulated_factor*str.fixed_displacement(:);str.phi - str.data.acumulated_factor*fixed_potential];
RT0_aux                                         =  RT0(str.freedof_prescribed);
RT0                                             =  RT0_aux;
else
%--------------------------------------------------------------------------
% RT0 according to rotational Dirichlet boundary conditions
%--------------------------------------------------------------------------
X_prescribed0                                   =  str.previous_Eulerian_x(str.freedof_prescribed);
X_prescribed0_                                  =  reshape(X_prescribed0,2,[]);
R                                               =  [[cos(theta)  -sin(theta)];[sin(theta) cos(theta)]];
RRT0                                            =  R*X_prescribed0_;
RRT0                                            =  reshape(RRT0,[],1);
RT0                                             =  str.Eulerian_x(str.freedof_prescribed) - RRT0;
end


Rx(str.freedof_prescribed)                      =  Rx(str.freedof_prescribed) + str.T0_arclength;   %  The residual Rx must include the contribution T0
Rx(str.fixdof_arclength,:)                      =  [];
RU                                              =  [Rx;RT0];
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Dx   =  -Ru - Rw0*Deltalambda 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Ru                                              =  (Kuu\RU);
Ru                                              =  Ru(1:size(Kxx,1));
 

total_size                                      =  size(physics.K_matrix,1);
Rw_size                                         =  size(w0_fixed_residual,1);
zeros_added                                     =  zeros(total_size-Rw_size,1);
%vector2                                        =  [zeros_added;w0_fixed_residual];
vector2                                         =  [p0_fixed_residual;w0_fixed_residual];
vector3                                         =  (Kuu\Kulambda);
vector3                                         =  vector3(1:size(Kxx,1));

Rw0                                             =  -Kxx\vector2 + vector3;


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Coefficients of the second order polynomial 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
freedof                                         =  str.freedof_arclength; 
switch str.arc_length_type
    case 'one_dof' 
         a                                      =  -Ru(constraint_dof_u);
         b                                      =  -Rw0(constraint_dof_u);
    case 'all_dof'
         am                                     =  -Ru(freedof<=str.n_nodes*3);
         bm                                     =  -Rw0(freedof<=str.n_nodes*3);
         ae                                     =  -Ru(freedof>str.n_nodes*3);
         be                                     =  -Rw0(freedof>str.n_nodes*3);
         a                                      =  [am*mechanical_scaling;ae*potential_scaling];
         b                                      =  [bm*mechanical_scaling;be*potential_scaling];
end
T                                               =  (total_u_fixed'*total_u_fixed + Rw0_fixed^2 + p0_fixed_residual'*p0_fixed_residual)*scaling^2;
A                                               =  b'*b + T;
B                                               =  2*a'*b + 2*solution.Du'*b + 2*solution.Dw0*(T);
C                                               =  a'*a + 2*a'*solution.Du + solution.Du'*solution.Du + (T)*solution.Dw0^2 - arc_length_radius^2;
Dw01                                            =  real(-B  -  sqrt(B^2 - 4*A*C))/(2*A);
Dw02                                            =  real(-B  +  sqrt(B^2 - 4*A*C))/(2*A);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
switch iteration 
    case 0
         %-----------------------------------------------------------------
         % Predictor according to criterion by Feng
         %-----------------------------------------------------------------
         if norm(solution.Du_old)>0
            Dw0                                 =  sign(solution.Du_old'*b)*arc_length_radius/sqrt(b'*b);
         else
            Dw01                                =  arc_length_radius/sqrt(b'*b);
            Dw02                                =  -arc_length_radius/sqrt(b'*b);
        
            Delta_solution1                     =  a + b*Dw01;
            Delta_solution2                     =  a + b*Dw02;
        
            switch str.arc_length_type
                   case 'one_dof'
                        u1                      =  Delta_solution1;
                        u2                      =  Delta_solution2;        
                   case 'all_dof'
                        u1                      =  Delta_solution1(constraint_dof_u);
                        u2                      =  Delta_solution2(constraint_dof_u);        
            end
            [maximum,identifier]                =  max([u1 u2]);
            switch identifier
                   case 1
                        Dw0                     =  Dw01;
                   case 2
                        Dw0                     =  Dw02;
            end                            
         end         
         %-----------------------------------------------------------------
         %  Update
         %-----------------------------------------------------------------
          solution.w0                           =  0;
          switch str.arc_length_type
               case 'all_dof'
                  solution.Du                   =  zeros(size(freedof,1),1);   
              case 'one_dof'
                  solution.Du                   =  0;
          end
          solution.Dw0                          =  0;
          solution.Dw0_incremental              =  Dw0;
    otherwise
          Delta_solution1                       =  a +  b*Dw01;  
          Delta_solution2                       =  a +  b*Dw02;
         

          sk1_1                                 =  [(solution.Du + Delta_solution1(1:size(solution.Du,1)));  (solution.Dw0 + Dw01)*Rw0_fixed*scaling];
          sk1_2                                 =  [(solution.Du + Delta_solution2(1:size(solution.Du,1)));  (solution.Dw0 + Dw02)*Rw0_fixed*scaling];
          cos_theta1                            =  (sk'*sk1_1)/arc_length_radius^2;
          cos_theta2                            =  (sk'*sk1_2)/arc_length_radius^2;
          theta1                                =  abs(acos(cos_theta1));
          theta2                                =  abs(acos(cos_theta2));
          [min_theta,identifier]                =  min([theta1 theta2]); 
         
          switch identifier 
                 case 1
                      solution.Dw0_incremental  =  Dw01;
                 case 2
                      solution.Dw0_incremental  =  Dw02; 
         end    
end     

%--------------------------------------------------------------------------
%  Update      
%--------------------------------------------------------------------------
solution.Du                                     =  solution.Du   +  a + b*solution.Dw0_incremental;
solution.DU                                     =  solution.Du;
solution.Dw0                                    =  solution.Dw0  +  solution.Dw0_incremental;
solution.w0                                     =  solution.w0   +  solution.Dw0_incremental;

%-----------------------------------------------------------------
%  Update of remaining variables (displacements, potential,
%  pressure... for the remaining nodels)    
%-----------------------------------------------------------------
str.solu.solution_incremented(str.freedof_arclength,...
    1)                                          =  -Ru - Rw0*solution.Dw0_incremental;
str                                             =  updated_incremental_potential(str);
str                                             =  updated_incremental_geometry(str);

%DT0                                             =  M*(RT0 - KxT0'*(Kxx\Rx) - X_prescribed0*solution.Dw0_incremental);
%DT0                                             =  -IT0T0*(Kuu\(RU + Kulambda*solution.Dw0_incremental));
DT0                                             =  -(Kuu\(RU + Kulambda*solution.Dw0_incremental));
DT0                                             =  DT0(size(Kxx,1)+1:size(DT0,1));
str.T0_arclength                                =  str.T0_arclength + DT0;


%--------------------------------------------------------------------------
%  Update acumulated factor    
%--------------------------------------------------------------------------
str.data.acumulated_factor                      =  str.data.acumulated_factor + solution.Dw0_incremental;
str.arc_length.theta                            =  theta + solution.Dw0_incremental;

%str                                             =  updated_constrained_potential(str);



