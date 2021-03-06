function [solution,str]                         =  output_arclength_v7(physics,solution,constraint_dof_u,Rw0_fixed,w0_fixed_residual,p0_fixed_residual,sk,total_u_fixed,str,iteration)


arc_length_radius                               =  str.arc_length_radius;

%mechanical_activation                          =  1;
%electric_activation                            =  1;
mechanical_scaling                              =  1;
%potential_scaling                              =  1e5*str.properties.material_parameters.e/sqrt(str.properties.material_parameters.mu*str.properties.material_parameters.e);
%potential_scaling                               =  str.properties.material_parameters.e(1)/sqrt(str.properties.material_parameters.mu(1)*str.properties.material_parameters.e(1));
potential_scaling                               =  1e1*str.properties.material_parameters.e(1)/sqrt(str.properties.material_parameters.mu(1)*str.properties.material_parameters.e(1));
%scaling                                        =  1*sqrt(str.properties.material_parameters.mu/str.properties.material_parameters.e)/str.data.scaling_area;

scaling                                         =  str.arc_length_scaling;

theta                                           =  str.arc_length.theta;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Stiffness matrices  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Kxx                                             =  physics.K_matrix;

%--------------------------------------------------------------------------
% Global assembly of the matrix
%--------------------------------------------------------------------------
Kuu                                             =  Kxx(str.freedof,str.freedof);
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
RU                                              =  Rx;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Dx   =  -Ru - Rw0*Deltalambda  
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
Ru                                              =  (Kuu\RU);
 

%vector2                                         =  [p0_fixed_residual;w0_fixed_residual];
elec_dofs                                       =  str.elec_freedof+str.n_nodes*str.data.dim;



Rauxiliar                                       =  zeros(size(Kxx,1),1);
Rauxiliar(elec_dofs,1)                          =  w0_fixed_residual;
Rauxiliar(str.fixdof,:)                         =  [];
Rauxiliar                                       =  -Kxx(str.freedof,str.freedof)\Rauxiliar;
Rw0                                             =  Rauxiliar;
%Rauxiliar                                       =  -Kxx(elec_dofs,elec_dofs)\w0_fixed_residual;                       




% Rw0                                             =  zeros(size(Kxx,1),1);
% Rw0(elec_dofs)                                  =  Rauxiliar;
% Rw0(str.fixdof)                                 =  [];

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


 
%--------------------------------------------------------------------------
%  Update acumulated factor     
%--------------------------------------------------------------------------
str.data.acumulated_factor                      =  str.data.acumulated_factor + solution.Dw0_incremental;
str.arc_length.theta                            =  theta + solution.Dw0_incremental;

%str                                             =  updated_constrained_potential(str);



