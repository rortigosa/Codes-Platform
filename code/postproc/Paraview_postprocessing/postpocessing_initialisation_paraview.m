function str                                  =  postpocessing_initialisation_paraview(str)

dim                                           =  str.data.dim;
nnode                                         =  str.postproc.n_nodes;
%--------------------------------------------------------------------------
% Eulerian coordinates and electric potential.
%--------------------------------------------------------------------------
str.postproc.Eulerian_x                       =  zeros(dim,nnode);
str.postproc.displacement                     =  zeros(dim,nnode);
str.postproc.phi                              =  zeros(nnode,1);
%--------------------------------------------------------------------------
% Cauchy stress and the different contributions.
%--------------------------------------------------------------------------
str.postproc.sigma                            =  zeros((dim+1)*dim/2,nnode);
str.postproc.sigma_pressure                   =  zeros(nnode,1);
%--------------------------------------------------------------------------
% SigmaF, SigmaH, SigmaJ.
%--------------------------------------------------------------------------
str.postproc.SigmaF                           =  zeros(dim,dim,nnode);
str.postproc.SigmaH                           =  zeros(dim,dim,nnode);
str.postproc.SigmaJ                           =  zeros(nnode,1);
str.postproc.First_Piola                      =  zeros(dim,dim,nnode);
%--------------------------------------------------------------------------
% Eulerian and Lagrangian electric field.
%--------------------------------------------------------------------------
str.postproc.E0                               =  zeros(dim,nnode);
str.postproc.E                                =  zeros(dim,nnode);
str.postproc.E_norm                           =  zeros(nnode,1);
str.postproc.electric_breakdown_factor        =  zeros(nnode,1);
%--------------------------------------------------------------------------
% Eulerian and Lagrangian electric displacement.
%--------------------------------------------------------------------------
str.postproc.D                                =  zeros(dim,nnode);
str.postproc.D0                               =  zeros(dim,nnode);
%--------------------------------------------------------------------------
% Deformation gradient, cofactor and jacobian.
%--------------------------------------------------------------------------
str.postproc.F                                =  zeros(dim,str.data.dim,nnode);
str.postproc.H                                =  zeros(dim,str.data.dim,nnode);
str.postproc.J                                =  zeros(nnode,1);
