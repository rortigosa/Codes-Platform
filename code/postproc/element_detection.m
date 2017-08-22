%**************************************************************************
% This function detects what element a node, given its coordinates, belongs
% to. (At the moment, it works correctly just for quadrilaterals.)
%**************************************************************************

function [element_sought]  =  element_detection(e_coordinates,str)

str.oldquadrature.Chi                     =  str.quadrature.Chi;
str.quadrature.Chi                        =  [0  0];
for ielem = 1:str.n_elem
    %---------------------------------------------
    % Needed variables
    %---------------------------------------------
    nodes_elem                            =  str.connectivity(ielem,:);
    xelem                                 =  str.Eulerian_x(:,nodes_elem);    
    Xelem                                 =  str.Lagrangian_X(:,nodes_elem);
    phielem                               =  zeros(str.n_node_elem,1);
    %---------------------------------------------
    % Center of the element in Lagrangian configuration.
    %---------------------------------------------
    center                                =  sum(str.Lagrangian_X(:,nodes_elem),2)/str.n_node_elem;
    %---------------------------------------------
    % Vector joining the center with the specified
    % coordinates in Lagrangian coordinates.
    %---------------------------------------------
    vector                                =  e_coordinates - center;
    %---------------------------------------------
    % Isoparametric coordinates of the vector.
    %---------------------------------------------
    [str]                                 =  gradients(xelem,Xelem,phielem,str);
    iso_vector                            =  inv(str.grad.DX_chi)*vector;
    
    iso_vector                            =  abs(iso_vector);
    if max(iso_vector)<1+1e-4
       element_sought                     =  ielem;
       break;
    end
end

str.quadrature.Chi                         =  str.oldquadrature.Chi;
