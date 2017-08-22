function [str]                                         =  connectivity_plot(str)

str.oldquadrature.Chi                                  =  str.quadrature.Chi;
str.quadrature.Chi                                     =  [0  0];  
for ielem=1:str.n_elem
    %---------------------------------------------
    % Needed variables
    %---------------------------------------------
    nodes_elem                                         =  str.connectivity(ielem,:);
    xelem                                              =  str.Eulerian_x(:,nodes_elem);    
    Xelem                                              =  str.Lagrangian_X(:,nodes_elem);
    phielem                                            =  zeros(str.n_node_elem,1);
    %---------------------------------------------
    % Center of the element in Lagrangian configuration.
    %---------------------------------------------
    center                                             =  sum(str.Lagrangian_X(:,nodes_elem),2)/str.n_node_elem;
    %---------------------------------------------
    % Vector joining the center with the specified
    % coordinates in Lagrangian coordinates.
    %---------------------------------------------
    vertex_counter                                     =  0;
    for inode=1:str.n_node_elem
        node                                           =  str.connectivity(ielem,inode);
        n_coord                                        =  str.Lagrangian_X(:,node);
        vector                                         =  n_coord - center;
        %---------------------------------------------
        % Isoparametric coordinates of the vector.
        %---------------------------------------------
        [str]                                          =  gradients(xelem,Xelem,phielem,str);
        iso_vector                                     =  inv(str.grad.DX_chi)*vector;
        %---------------------------------------------
        % Detecting if the node is a vertex node or not.
        %---------------------------------------------
        [vertex_node]                                  =  vertex_nodes_detection(str,iso_vector);
        if vertex_node==1
           vertex_counter                              =  vertex_counter + 1;
           str.connectivity_plot(ielem,vertex_counter) =  node;
           str.localconnect_plot(ielem,vertex_counter) =  inode;
        end
    end
    str.auxiliar                                       =  1:1:str.n_node_elem;
    v                                                  =  str.localconnect_plot(ielem,:);
    str.auxiliar(v)                                    =  [];
    str.deleted_localnode(ielem,:)                     =  str.auxiliar;    
end


str.quadrature.Chi                                      =  str.oldquadrature.Chi;

end