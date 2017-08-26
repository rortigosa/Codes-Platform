function  ParaviewPostprocessor(str,filename)

local_ordering                        =  local_ordering_nodes(2);                                                            
n_elem                                =  str.mesh.volume.n_elem;
n_nodes                               =  (8*1)*str.fem.postprocessing.degree^3*n_elem;
connectivity                          =  reshape((1:n_nodes)',8,[]);
%--------------------------------------------------------------------------
%  Initialise postprocessing variables 
%--------------------------------------------------------------------------
x_solid                               =  zeros(3,n_nodes);
X_solid                               =  zeros(3,n_nodes);
F_solid                               =  zeros(9,n_nodes);
H_solid                               =  zeros(9,n_nodes);
J_solid                               =  zeros(n_nodes,1);
p_solid                               =  zeros(n_nodes,1);
Piola_solid                           =  zeros(9,n_nodes);
magnifying_factor                     =  1e2;
str.solution.x0.Eulerian_x            =  str.solution.x0.Eulerian_x;
str.solution.x0.Lagrangian_X          =  str.solution.x0.Lagrangian_X;
%--------------------------------------------------------------------------
%  Loop over elements 
%--------------------------------------------------------------------------
fem.x0                                =  str.fem.x0.postprocessing;
fem.pressure                          =  str.fem.pressure.postprocessing;
%new_elem                             =  1;
final                                 =  0;
for ielem=1:n_elem       
    n_gauss_points                    =  str.mesh.x0.n_node_elem;
    %----------------------------------------------------------------------
    % x and X of the postprocessing mesh mesh 
    %----------------------------------------------------------------------
    xelem                             =  zeros(3,size(str.fem.x0.nodes.N,2));
    Xelem                             =  zeros(3,size(str.fem.x0.nodes.N,2));
    for inode=1:size(str.fem.x0.nodes.N,2)
        xelem(:,inode)                =  str.solution.x0.Eulerian_x(:,str.mesh.x0.connectivity(:,ielem))*fem.x0.N(:,inode);
        Xelem(:,inode)                =  str.solution.x0.Lagrangian_X(:,str.mesh.x0.connectivity(:,ielem))*fem.x0.N(:,inode);        
    end
    %----------------------------------------------------------------------
    % Gradients in the nodes of the postprocessing mesh
    %----------------------------------------------------------------------
    gradients                         =  gradients_shell(ielem,n_gauss_points,fem,str.mesh,str.solution);
    %----------------------------------------------------------------------
    % Obtain pressure at every node of the postprocessing mesh
    %----------------------------------------------------------------------
    pressure                          =  pressure_gauss_point(n_gauss_points,fem.pressure,str.solution.pressure(str.mesh.pressure.connectivity(:,ielem)));
    %----------------------------------------------------------------------
    % Pre-stretch at every node of the postprocessing mesh
    %----------------------------------------------------------------------
    Pre_stretch                       =  Pre_stretch_gauss_point_final(str.material_information.layer(ielem),n_gauss_points,str.material_information.material_parameters,gradients.x0.Normal_vector,str.time_integrator.time_iteration);
    %----------------------------------------------------------------------
    % First and second derivatives of the model  at every node of the
    % postprocdessing mesh
    %----------------------------------------------------------------------
    DU                                =  first_derivative_computation(gradients.F,gradients.H,gradients.J,...
                                                                    n_gauss_points,str.material_information.material_model,str.material_information.material_identifier(ielem),...
                                                                    str.material_information.material_parameters,Pre_stretch);
    %----------------------------------------------------------------------
    % First Piola-Kirchhoff stress tensor.
    %----------------------------------------------------------------------
    Piola                             =  First_Piola_Kirchhoff_stress_tensor_final(n_gauss_points,gradients,pressure,DU);
    %----------------------------------------------------------------------
    % First Piola-Kirchhoff stress tensor.  
    %----------------------------------------------------------------------
    for idegree=1:str.fem.postprocessing.degree^3
        %new_elem                     =  new_elem + 1;
        initial                       =  final + 1;
        final                         =  final + 8;
 
        %connectivity(:,new_elem)     =  
        x_solid(:,initial:final)      =  magnifying_factor*xelem(:,local_ordering(:,idegree));
        X_solid(:,initial:final)      =  magnifying_factor*Xelem(:,local_ordering(:,idegree));
        p_solid(initial:final)        =  pressure(local_ordering(:,idegree));
        F_solid(:,initial:final)      =  reshape(gradients.F(:,:,local_ordering(:,idegree)),9,[]);
        H_solid(:,initial:final)      =  reshape(gradients.H(:,:,local_ordering(:,idegree)),9,[]);
        J_solid(initial:final)        =  gradients.J(local_ordering(:,idegree));
        Piola_solid(:,initial:final)  =  reshape(Piola(:,:,local_ordering(:,idegree)),9,[]);
    end
end

ux                                    =  reshape(x_solid(1,:) - X_solid(1,:),[],1);
uy                                    =  reshape(x_solid(2,:) - X_solid(2,:),[],1);
uz                                    =  reshape(x_solid(3,:) - X_solid(3,:),[],1);
Fxx                                   =  reshape(F_solid(1,:),[],1);
Fyx                                   =  reshape(F_solid(2,:),[],1);
Fzx                                   =  reshape(F_solid(3,:),[],1);
Fxy                                   =  reshape(F_solid(4,:),[],1);
Fyy                                   =  reshape(F_solid(5,:),[],1);
Fzy                                   =  reshape(F_solid(6,:),[],1);
Fxz                                   =  reshape(F_solid(7,:),[],1);
Fyz                                   =  reshape(F_solid(8,:),[],1);
Fzz                                   =  reshape(F_solid(9,:),[],1);
Hxx                                   =  reshape(H_solid(1,:),[],1);
Hyx                                   =  reshape(H_solid(2,:),[],1);
Hzx                                   =  reshape(H_solid(3,:),[],1);
Hxy                                   =  reshape(H_solid(4,:),[],1);
Hyy                                   =  reshape(H_solid(5,:),[],1);
Hzy                                   =  reshape(H_solid(6,:),[],1);
Hxz                                   =  reshape(H_solid(7,:),[],1);
Hyz                                   =  reshape(H_solid(8,:),[],1);
Hzz                                   =  reshape(H_solid(9,:),[],1);
J                                     =  reshape(J_solid,[],1);
p                                     =  reshape(p_solid,[],1);
Pxx                                   =  reshape(Piola_solid(1,:),[],1);
Pyx                                   =  reshape(Piola_solid(2,:),[],1);
Pzx                                   =  reshape(Piola_solid(3,:),[],1);
Pxy                                   =  reshape(Piola_solid(4,:),[],1);
Pyy                                   =  reshape(Piola_solid(5,:),[],1);
Pzy                                   =  reshape(Piola_solid(6,:),[],1);
Pxz                                   =  reshape(Piola_solid(7,:),[],1);
Pyz                                   =  reshape(Piola_solid(8,:),[],1);
Pzz                                   =  reshape(Piola_solid(9,:),[],1);
  
VTK_plotting(filename,'unstructured_grid',x_solid',connectivity','scalars','ux',ux,...
                                                                'scalars','uy',uy,...
                                                                'scalars','uz',uz,...
                                                                'scalars','Fxx',Fxx,...
                                                                'scalars','Fxy',Fxy,...
                                                                'scalars','Fxz',Fxz,...
                                                                'scalars','Fyx',Fyx,...
                                                                'scalars','Fyy',Fyy,...
                                                                'scalars','Fyz',Fyz,...
                                                                'scalars','Fzx',Fzx,...
                                                                'scalars','Fzy',Fzy,...
                                                                'scalars','Fzz',Fzz,...
                                                                'scalars','Hxx',Hxx,...
                                                                'scalars','Hxy',Hxy,...
                                                                'scalars','Hxz',Hxz,...
                                                                'scalars','Hyx',Hyx,...
                                                                'scalars','Hyy',Hyy,...
                                                                'scalars','Hyz',Hyz,...
                                                                'scalars','Hzx',Hzx,...
                                                                'scalars','Hzy',Hzy,...
                                                                'scalars','Hzz',Hzz,...
                                                                'scalars','J',J,...
                                                                'scalars','p',p,...
                                                                'scalars','Pxx',Pxx,...
                                                                'scalars','Pxy',Pxy,...
                                                                'scalars','Pxz',Pxz,...
                                                                'scalars','Pyx',Pyx,...
                                                                'scalars','Pyy',Pyy,...
                                                                'scalars','Pyz',Pyz,...
                                                                'scalars','Pzx',Pzx,...
                                                                'scalars','Pzy',Pzy,...
                                                                'scalars','Pzz',Pzz) 

%end                                                            
                                                            
function local_ordering               =  local_ordering_nodes(degree)                                                            
ordering1                             =  zeros(8,1);
inode                                 =  1;
for iz=1:2
    for iy=1:2
        for ix=1:2
            ordering1(inode)          =  ix + (iy - 1)*(degree + 1) + (iz - 1)*(degree + 1)^2;
            inode                     =  inode + 1;
        end
    end
end
ordering2                             =  zeros(degree^3,1);
inode                                 =  1;
for iz=1:degree
    for iy=1:degree
        for ix=1:degree
            ordering2(inode)          =  ix + (iy - 1)*(degree + 1) + (iz - 1)*(degree + 1)^2;
            inode                     =  inode + 1;
        end
    end
end
local_ordering                        =  zeros(8,str.fem.x0.degree);
for ielem=1:degree^3           
    local_ordering(:,ielem)           =  ordering1 + (ordering2(ielem) - 1);
end
local_ordering                        =  local_ordering([1 2 4 3 5 6 8 7],:);    
end
end