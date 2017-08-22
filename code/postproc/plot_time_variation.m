function    [str]  = plot_time_variation(str)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%  Beam formulation
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% switch str.data.dim
%     case {13,23}
% X0                                                  =  [0.5;0.5;6];
% X0_discrete                                         =  [0;0;6];
% X0_center                                           =  [0;0;5];
% for inode=1:size(str.solid.BEAM_SHELL.continuum.Eulerian_x,2)
%     X                                               =  str.solid.BEAM_SHELL.continuum.Lagrangian_X(:,inode);
%     r                                               =  X0-X;
%     if norm(r)<1e-7
%        node                                         =  inode;
%        break;
%     end
% end
% 
% for inode=1:size(str.solid.BEAM_SHELL.discrete.Eulerian_x,2)
%     X_discrete                                      =  str.solid.BEAM_SHELL.discrete.Lagrangian_X(:,inode);
%     r                                               =  X0_discrete-X_discrete;
%     if norm(r)<1e-7
%        node_discrete                                =  inode;
%        break;
%     end
% end
% 
% for inode=1:size(str.solid.BEAM_SHELL.discrete.Eulerian_x,2)
%     X_center                                        =  str.solid.BEAM_SHELL.discrete.Lagrangian_X(:,inode);
%     r                                               =  X0_center-X_center;
%     if norm(r)<1e-7
%        center                                =  inode;
%        break;
%     end
% end
% 
% figure(3)
% str.post.displ1(str.time_iteration)                  =  norm(str.solid.BEAM_SHELL.continuum.Eulerian_x(:,node) - ...
%                                                              str.solid.BEAM_SHELL.continuum.Lagrangian_X(:,node));
% 
% str.post.displ2(str.time_iteration)                  =  (str.solid.BEAM_SHELL.continuum.Eulerian_x(2,node) - ...
%                                                              str.solid.BEAM_SHELL.continuum.Lagrangian_X(2,node));
% 
% str.post.displ3(str.time_iteration)                  =  (str.solid.BEAM_SHELL.continuum.Eulerian_x(3,node) - ...
%                                                              str.solid.BEAM_SHELL.continuum.Lagrangian_X(3,node));
% %str.post.displ1(str.time_iteration)                  =  str.solid.BEAM_SHELL.discrete.Eulerian_x(1,node_discrete) - ...
% %                                                             str.solid.BEAM_SHELL.discrete.Lagrangian_X(1,node_discrete);
% 
% %str.post.displ2(str.time_iteration)                  =  str.solid.BEAM_SHELL.discrete.Eulerian_x(2,node_discrete) - ...
% %                                                             str.solid.BEAM_SHELL.discrete.Lagrangian_X(2,node_discrete);
% % str.post.displ2(str.time_iteration)                  =  str.solid.BEAM_SHELL.discrete.Eulerian_x(1,center) - ...
% %                                                              str.solid.BEAM_SHELL.discrete.Lagrangian_X(1,center);
%                                                          
%                                                          
%                                                          v=reshape(str.solid.BEAM_SHELL.discrete.velocity,6,size(str.solid.BEAM_SHELL.discrete.velocity,1)/6);
% a=reshape(str.solid.BEAM_SHELL.discrete.acceleration,6,size(str.solid.BEAM_SHELL.discrete.acceleration,1)/6);
% 
% linear_v  =  v(1:3,:);
% angular_v =  v(4:6,:);
% linear_a  =  a(1:3,:);
% angular_a =  a(4:6,:);
% 
% X_continuum =  str.solid.BEAM_SHELL.continuum.Eulerian_x;
% 
% % str.post.displ2(str.time_iteration)                  =  norm(linear_v(:,node_discrete));
% % str.post.displ3(str.time_iteration)                  =  norm(linear_a(:,node_discrete));
% % str.post.displ4(str.time_iteration)                  =  norm(angular_v(:,node_discrete));
% % str.post.displ5(str.time_iteration)                  =  norm(angular_a(:,node_discrete));
% 
% 
% str.post.displ4(str.time_iteration)                  =  linear_v(3,node_discrete);
% str.post.displ5(str.time_iteration)                  =  linear_a(3,node_discrete);
% str.post.displ6(str.time_iteration)                  =  angular_v(3,node_discrete);
% str.post.displ7(str.time_iteration)                  =  angular_a(3,node_discrete);
% 
% 
% subplot(3,3,1);  
% plot(str.post.displ1(1:str.time_iteration),'-*')
% title('discrete (beam) displacement at particular point')
% 
% 
% subplot(3,3,2);  
% plot(str.post.displ2(1:str.time_iteration),'-*')
% title('y (beam) displacement center of the beam')
% 
% subplot(3,3,3);  
% plot(str.post.displ3(1:str.time_iteration),'-*')
% title('z (beam) displacement center of the beam')
% 
% 
% subplot(3,3,4); 
% plot(str.post.displ4(1:str.time_iteration),'-*')
% title('norm of discrete (beam) linear velocity at particular point')
% 
% subplot(3,3,5); 
% plot(str.post.displ5(1:str.time_iteration),'-*')
% title('norm of discrete (beam) linear acceleration at particular point')
% 
% subplot(3,3,6); 
% plot(str.post.displ6(1:str.time_iteration),'-*')
% title('norm of discrete (beam) angular velocity at particular point')
% 
% subplot(3,3,7); 
% plot(str.post.displ7(1:str.time_iteration),'-*')
% title('norm of discrete (beam) angular acceleration at particular point')
% 
% subplot(3,3,8); 
% plot3(X_continuum(1,:),X_continuum(2,:),X_continuum(3,:),'o')
% title('Deformed configuration')
% 
% 
% end
% 
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %  Continuum formulation
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% switch str.data.dim
%     case {2,3}
% X0                                                  =  [0.5;0.5;6];
% for inode=1:str.n_nodes
%     X                                               =  str.nodes(inode,:)';
%     r                                               =  X0-X;
%     if norm(r)<1e-7
%        node                                         =  inode;
%        break;
%     end
% end
% 
%--------------------------------------------------------------------------
% Energy conservation
%--------------------------------------------------------------------------
[Kinetic_energy,Elastic_energy,Total_energy,...
    Kinetic_power,Elastic_power,Total_power]         =  conservation_energy_monitoring(str);

figure(3)
% str.post.displ1(str.time_iteration)                  =  norm(str.Eulerian_x(:,node) - str.Lagrangian_X(:,node));
% v=reshape(str.velocity,3,size(str.velocity,1)/3);
% a=reshape(str.acceleration,3,size(str.acceleration,1)/3);
% str.post.displ2(str.time_iteration)                  =  norm(v(:,node));
% str.post.displ3(str.time_iteration)                  =  (a(2,node));
str.post.Kinetic_energy(str.time_iteration)          =  Kinetic_energy;
str.post.Elastic_energy(str.time_iteration)          =  Elastic_energy;
str.post.Total_energy(str.time_iteration)            =  Total_energy;
str.post.Kinetic_power(str.time_iteration)           =  Kinetic_power;
str.post.Elastic_power(str.time_iteration)           =  Elastic_power;
str.post.Total_power(str.time_iteration)             =  Total_power;

% subplot(3,3,1);  
% plot(str.post.displ1(1:str.time_iteration),'-*')
% title('displacement at particular point')
% 
% subplot(3,3,2); 
% plot(str.post.displ2(1:str.time_iteration),'-*')
% title('norm of velocity at particular point')
% 
% subplot(3,3,3); 
% plot(str.post.displ3(1:str.time_iteration),'-*')
% title('norm of acceleration at particular point')

subplot(3,3,4); 
plot(str.post.Kinetic_energy(1:str.time_iteration),'-*')
title('Kinetic energy')

subplot(3,3,5); 
plot(str.post.Elastic_energy(1:str.time_iteration),'-*')
title('Elastic energy')

subplot(3,3,6); 
plot(str.post.Total_energy(1:str.time_iteration),'-*')
title('Total energy')

subplot(3,3,7); 
plot(str.post.Kinetic_power(1:str.time_iteration),'-*')
title('Kinetic power')

subplot(3,3,8); 
plot(str.post.Elastic_power(1:str.time_iteration),'-*')
title('Elastic power')

subplot(3,3,9); 
plot(str.post.Total_power(1:str.time_iteration),'-*')
title('Total power')

end
% 
% 
% 
% %str.plotting.phi_node                              =  241;
% % str.plotting.phi_node                               =  66;
% % str.post.phi(1,1)                                   =  0;
% % str.post.phi(str.time_iteration)                    =  str.phi(str.plotting.phi_node,1);
% % figure(200)
% % plot(str.post.phi(1:str.time_iteration),'-*')
