%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes:
% A. The interpolated values at every Gauss point for the unknowns of the
% particular mixed formulation.
% B. With the information in A, in calculates the partial derivatives of
% the particular stored energy functional associated to the considered
% mixed formulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   gauss_level_information                        =  postprocessing_element_nodes_paraview(str)


ielem                                                     =  str.ielem;
gauss_level_information                                   =  [];

nodes_x                                                   =  str.connectivity(ielem,:);
%nodes_x                                                  =  str.postproc.connectivity(ielem,:);

nodes_phi                                                 =  str.connectivity(ielem,:);
%nodes_phi                                                =  str.postproc.connectivity(ielem,:);
nodal_level_information.Lagrangian_X                      =  str.Lagrangian_X(:,nodes_x);
nodal_level_information.Eulerian_x                        =  str.Eulerian_x(:,nodes_x);

velocity                                                  =  reshape(str.velocity,str.data.dim,str.n_nodes);
nodal_level_information.velocity                          =  velocity(:,nodes_x);
acceleration                                              =  reshape(str.acceleration,str.data.dim,str.n_nodes);
nodal_level_information.acceleration                      =  acceleration(:,nodes_x);
nodal_level_information.phi                               =  str.phi(nodes_phi);

for igauss=1:size(str.quadrature.Chi,1)
    %----------------------------------------------------------------------
    % Shape functions
    %----------------------------------------------------------------------
    Nu                                                    =  str.f_e.N(:,igauss);
    Nphi                                                  =  str.f_e.N(:,igauss);
    %----------------------------------------------------------------------
    % Information at Gauss level
    %----------------------------------------------------------------------
    gauss_level_information.Lagrangian_X(:,igauss)        =  nodal_level_information.Lagrangian_X*Nu;    
    gauss_level_information.Eulerian_x(:,igauss)          =  nodal_level_information.Eulerian_x*Nu;    
    gauss_level_information.velocity(:,igauss)            =  nodal_level_information.velocity*Nu;
    gauss_level_information.acceleration(:,igauss)        =  nodal_level_information.acceleration*Nu;
    
    
    gauss_level_information.phi(igauss)                   =  nodal_level_information.phi'*Nphi;
end
