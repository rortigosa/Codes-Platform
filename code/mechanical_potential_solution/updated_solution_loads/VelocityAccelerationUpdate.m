
function solution                         =  VelocityAccelerationUpdate(time_integrator,solution)

switch time_integrator.type
    case {'Newmark_beta','generalised_alpha'}
         beta                             =  time_integrator.beta;
         gamma                            =  time_integrator.gamma;
         Dt                               =  time_integrator.Deltat;
         v_factor                         =  gamma/(beta*Dt);
         a_factor                         =  1/(beta*Dt^2);
         
         solution.x.Eulerian_x_dot       =  solution.x.Eulerian_x_dot + D_x*v_factor;    
         solution.x.Eulerian_x_dot_dot   =  solution.x.Eulerian_x_dot_dot + D_x*a_factor;    
end
