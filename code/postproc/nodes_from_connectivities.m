function [str]                               =  nodes_from_connectivities(str)
nnodes                                       =  str.n_nodes;
str.n_nodes                                  =  0;
connectivity_vector                          =  reshape(str.connectivity',1,size(str.connectivity,2)*size(str.connectivity,1));

for inode=1:nnodes
      for icounter = 1:size(connectivity_vector,2)
          if connectivity_vector(icounter) == inode
             str.n_nodes                     =  str.n_nodes + 1;
             elem                            =  ceil(icounter/size(str.connectivity,2));
             local_node                      =  icounter - size(str.connectivity,2)*(elem-1);
             node                            =  str.connectivity(elem,local_node);
             str.new_nodes(str.n_nodes,:)    =  str.nodes(node,:); 
             str.new_node_numb(str.n_nodes)  =  node;
             break;
          end
      end     
end
 
str.newdeleted_node_numb                     =  1:nnodes;
str.newdeleted_node_numb(str.new_node_numb)  =  [];

