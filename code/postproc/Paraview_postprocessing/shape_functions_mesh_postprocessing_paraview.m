%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  This function re-computes connectivities and shape functions due to the
%  inconsistency between the connectivities of both computationand and 
%  postprocessing meshes.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                           =  shape_functions_mesh_postprocessing_paraview(str)


str.f_e                                =  str.solution_nodes.f_e;
str.quadrature.Chi                     =  str.solution_nodes.quadrature.Chi;
% switch str.data.shape 
%     case 1
%          renumbering                       =  [1 2 4 3 5 6 8 7];
%          str.f_e  =  str.solution_nodes.f_e.N(renumbering,renumbering);
%          str.quadrature.Chi                =  str.solution_nodes.quadrature.Chi(renumbering,:);
%          str.postproc.connectivity     =  str.postproc.connectivity(:,renumbering);
%          str.connectivity              =  str.connectivity(:,renumbering);
%            %--------------
%            % Old postproc
%            %--------------
%          %renumbering                       =  [1 2 4 3 5 6 8 7];
% %          str.f_e  =  str.solution_nodes.f_e.N;
% %          str.quadrature.Chi                =  str.solution_nodes.quadrature.Chi;
% %          str.postproc.connectivity     =  str.postproc.connectivity;
% %          str.connectivity              =  str.connectivity;
% 
% 
% end
         
str.nodes_counter                                =  zeros(1,size(str.postproc.nodes,1));
for ielem=1:size(str.postproc.connectivity,1)
    for inode=1:size(str.postproc.connectivity,2)        
        node                           =  str.postproc.connectivity(ielem,inode);
        str.nodes_counter(1,node)      =  str.nodes_counter(1,node) + 1;
    end
end
    
  
% for ielem=1:str.n_elem
%     for inode=1:str.postproc.n_node_elem
%         node                                     =  str.postproc.connectivity(ielem,inode);
%         str.nodes_counter(1,node)                =  str.nodes_counter(1,node) + 1;
%     end
% end
  