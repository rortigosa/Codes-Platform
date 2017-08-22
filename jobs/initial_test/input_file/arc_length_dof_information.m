function [coordinate,constrained_dofs,...
    scaling_area,str]              =   arc_length_dof_information(str)

coordinate                        =  [0;  0;   10];
constrained_dofs                  =  [0   1];  %  dof 3 for displacements and dof 1 (redundant) for electric charge

 value                             =  [3e-3;3e-3];
 loc                               =  [0;0.05];      %  The surface charge is in coordinate = 0.5 (where coordinate is aligned with the cooridinate associated to degree) 
 degree                            =  [3;3];
 scaling_area                      =  36;




switch str.data.shape
    case 0
         n_nodes_face                                =  sum(1:(str.data.degree+1));
    case 1
         n_nodes_face                                =  (str.data.degree+1)^2;
end

elem                                                 =  1;
for ielem=1:str.n_elem
    for icase=1:length(value)
        element_nodes                                =  str.connectivity(ielem,:);    
        element_coordinates                          =  abs(str.nodes(element_nodes,degree(icase)) - loc(icase));    
    
       [element_coordinates,order]                   =  sort(element_coordinates);
       coordinates_to_consider                       =  element_coordinates(1:n_nodes_face);
       criterion                                     =  norm(coordinates_to_consider);
       %-------------------------------------------------------------------
       % X_volume
       %-------------------------------------------------------------------
       if max(str.nodes(element_nodes,1))-5<1e-6 
          if max(str.nodes(element_nodes,degree(icase)))-0.025<1e-6 
            if norm(criterion)<1e-6
               local_nodes                                 =  element_nodes(order(1:n_nodes_face));
               str.solid.bc.surface_charge.element(elem,:) =  [ielem  value(icase) local_nodes];
               elem                                        =  elem + 1;
            end
          end
       end
       if min(str.nodes(element_nodes,1))-5>1e-6 
          if max(str.nodes(element_nodes,degree(icase)))-0.025>1e-6 
            if norm(criterion)<1e-6
               local_nodes                                 =  element_nodes(order(1:n_nodes_face));
               str.solid.bc.surface_charge.element(elem,:) =  [ielem  value(icase) local_nodes];
               elem                                        =  elem + 1;
            end
          end
       end
    end
end
