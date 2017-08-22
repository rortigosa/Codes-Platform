function  BEAM_SHELL_MATLAB_plottings(str)

for imechanical_node=1:size(str.solid.BEAM_SHELL.discrete.Lagrangian_X,2)
    X_coor                    =  str.solid.BEAM_SHELL.discrete.Eulerian_x(:,imechanical_node);
    %information              =   imechanical_node;
    plot3(X_coor(1),X_coor(2),X_coor(3),'ro')
    %text_plot_information(X_coor,information)
    hold on
end
   
%for inode=49:size(str.solid.BEAM_SHELL.continuum.nodes)

for inode=1:size(str.solid.BEAM_SHELL.continuum.nodes)
    X_coor                    =  str.solid.BEAM_SHELL.continuum.Eulerian_x(:,inode);
    %information              =  inode;
    plot3(X_coor(1),X_coor(2),X_coor(3),'bo')
    %text_plot_information(X_coor,information)
    hold on
end


for imechanical_element=1:size(str.solid.BEAM_SHELL.discrete.mesh.connectivities,1)
    for imechanical_node=1:size(str.solid.BEAM_SHELL.discrete.mesh.connectivities,2)
        triad                 =  str.solid.BEAM_SHELL.discrete.Eulerian_covariant(:,:,imechanical_node,imechanical_element);
        for itriad=1:3
            vector            =  triad(:,itriad);
            mechanical_node   =  str.solid.BEAM_SHELL.discrete.mesh.connectivities(imechanical_element,imechanical_node);
            position          =  str.solid.BEAM_SHELL.discrete.Eulerian_x(:,mechanical_node);
            vector_plot(position,vector)
            hold on
        end
    end
end

            
        
    
    