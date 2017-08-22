function    [str]                                                = plot_arclength_variation(str)
%str.post.displ(str.acr_length.arc_length_increment)              =  str.Lagrangian_X(2,12)-str.Eulerian_x(2,12);
%str.post.displ(str.acr_length.arc_length_increment)              =  str.Lagrangian_X(2,2)-str.Eulerian_x(2,2);
str.post.displ(str.arc_length.arc_length_increment)              =  str.Lagrangian_X(2,187)-str.Eulerian_x(2,187);
str.post.displ(str.arc_length.arc_length_increment)              =  str.Lagrangian_X(2,8)-str.Eulerian_x(2,8);
figure(1)

%str.plot.yvariable(str.acr_length.arc_length_increment)          =  abs(str.nodal_loads.P(24,1))*str.arc_length.lambda;
%str.plot.yvariable(str.acr_length.arc_length_increment)          =  abs(str.nodal_loads.P(4,1))*str.arc_length.lambda;
str.plot.yvariable(str.arc_length.arc_length_increment)          =  abs(str.phi(2,1));
str.plot.xvariable                                               =  str.post.displ;
%as=str.plot.xvariable
plot(str.plot.xvariable,str.plot.yvariable,'-b')
