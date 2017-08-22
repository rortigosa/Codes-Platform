function str                                  =  postpocessing_acoustic_tensor_initialisation(str)

dim                                           =  str.data.dim;
nnode                                         =  str.postproc.n_nodes;
%-----------------------------------------
% Eulerian coordinates and electric potential.
%-----------------------------------------
str.postproc.Eulerian_x                       =  zeros(dim,nnode);
%-----------------------------------------
% Second Piola and the different contributions.
%-----------------------------------------
str.postproc.q                                =  zeros(nnode,1);
