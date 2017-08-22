%**************************************************************************
% This function 
%**************************************************************************

function [str]  =  plot_variables_treatment(str)

%---------------------------------------------------------
% New connectivities in the submesh.
%---------------------------------------------------------
%str.connectivity                          =  str.connectivity_plot;
%str.postproc.t                            =  str.postproc.post_connectivity';
str.postproc.t                            =  str.postproc.connectivity';

%---------------------------------------------------------
% Identify the nodes of the submesh and delete the ones the
% previous mesh
%---------------------------------------------------------
%[str]                                     =  nodes_from_connectivities(str);
%str.nodes                                 =  str.new_nodes;
str.postproc.nodes                        =  str.postproc.nodes';

%---------------------------------------------------------
% Removed the value of the variables in the mid nodes.
%---------------------------------------------------------
% v                                         =  str.newdeleted_node_numb';
% str.Lagrangian_X(:,v)                     =  [];
% str.Eulerian_x(:,v)                       =  [];
% str.phi(v)                                =  [];
% str.postproc.stress(:,v)                  =  [];
% str.postproc.sigma(:,v)                   =  [];
% str.postproc.D(:,v)                       =  [];
% str.postproc.D0(:,v)                      =  [];
% str.postproc.E0(:,v)                      =  [];
% str.postproc.E(:,v)                       =  [];
% str.postproc.F(:,:,v)                     =  [];
% str.nodes_counter(v)                      =  [];
% 
% v_vector                                  =  zeros(str.data.dim*length(v),1);
% compv                                     =  length(v);
% for iloop=1:str.data.dim
%     final                                 =  str.data.dim*compv - (str.data.dim - iloop);  
%     v_vector(iloop:str.data.dim:final,1)  =  str.data.dim*v - (str.data.dim - iloop)*ones(length(v),1);
% end
% str.assemb_force.total_force(v_vector)    =  [];
% 
% 



