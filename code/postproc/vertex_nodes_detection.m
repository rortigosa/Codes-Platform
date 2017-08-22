function [vertex_node]        =  vertex_nodes_detection(str,iso_vector)

vertex_node                   =  0;

switch str.data.dim
    case 2
         if norm(iso_vector - [-1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [-1;1])<1e-3
             vertex_node                =  1;
         end
    case 3
         if     norm(iso_vector - [-1;-1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;-1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [-1;1;-1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [-1;-1;1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;-1;1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [1;1;1])<1e-3
             vertex_node                =  1;
         elseif norm(iso_vector - [-1;1;1])<1e-3
             vertex_node                =  1;
         end
end        
