%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  In this function we compute different variables for a specific point.
%  Valid for mixed formulations. Not very generic. We net to load first the
%  mat file and then we can run it.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function variables  =  postprocessing_variables_gauss_node(X_gauss,str)

structure                                     =  str;
structure.data.degree                         =  str.data.degree;
structure                                     =  gauss_quadrature(structure);
structure.f_e.N                               =  [];
structure.f_e.DN_chi                          =  [];
structure.f_e.isoparametric_nodal_positions   =  [];
structure_u                                   =  shape_function_displacement(structure);

structure                                     =  str;
structure.data.degree                         =  str.data.degree_F;
structure                                     =  gauss_quadrature(structure);
structure.f_e.N                               =  [];
structure.f_e.DN_chi                          =  [];
structure.f_e.isoparametric_nodal_positions   =  [];
structure_F                                   =  shape_function_displacement(structure);

structure                                     =  str;
structure.data.degree                         =  str.data.degree_H;
structure                                     =  gauss_quadrature(structure);
structure.f_e.N                               =  [];
structure.f_e.DN_chi                          =  [];
structure.f_e.isoparametric_nodal_positions   =  [];
structure_H                                   =  shape_function_displacement(structure);

structure                                     =  str;
structure.data.degree                         =  str.data.degree_J;
structure                                     =  gauss_quadrature(structure);
structure.f_e.N                               =  [];
structure.f_e.DN_chi                          =  [];
structure.f_e.isoparametric_nodal_positions   =  [];
structure_J                                   =  shape_function_displacement(structure);

str.f_e.symbolic_N_u                          =  structure_u.f_e.symbN;
str.f_e.symbolic_DN_chi_u                     =  structure_u.f_e.symbDN_chi;
str.f_e.symbolic_N_F                          =  structure_F.f_e.symbN;
str.f_e.symbolic_DN_chi_F                     =  structure_F.f_e.symbDN_chi;
str.f_e.symbolic_N_H                          =  structure_H.f_e.symbN;
str.f_e.symbolic_DN_chi_H                     =  structure_H.f_e.symbDN_chi;
str.f_e.symbolic_N_J                          =  structure_J.f_e.symbN;
str.f_e.symbolic_DN_chi_J                     =  structure_J.f_e.symbDN_chi;

for ielem_coarse=1:str.n_elem
    %--------------------------------------------------------------
    % Try to detect what element of the coarse mesh
    % the Gauss point belong to.
    %--------------------------------------------------------------
    nodes                                     =  str.connectivity(ielem_coarse,:);
    Xelem                                     =  str.Lagrangian_X(:,nodes);
    %--------------------------------------------------------------
    % We estimate if the element contains the node.
    %--------------------------------------------------------------
    Characteristic_Cube                       =  str.mesh_information.Characteristic_Cube(:,:,ielem_coarse) ;
    [estimated_belonging]                     =  estimated_belonging_function(X_gauss,Characteristic_Cube,1.1);
    if estimated_belonging
        [iso_coor]                            =  Lagrangian_node_iso_domain(str,Xelem,X_gauss);
        [element_belonging]                   =  node_element_belonging(iso_coor,str);
        %-----------------------------------------------------------
        % Detection case. 
        %-----------------------------------------------------------
        if element_belonging == 1
            %--------------------------------------------------------
            % Evaluation of shape functions at the particular gauss node of the coarse mesh.
            %--------------------------------------------------------
            nodes                             =  str.connectivity(ielem_coarse,:);
            Xelem                             =  str.Lagrangian_X(:,nodes);
            xelem                             =  str.Eulerian_x(:,nodes);
            phielem                           =  str.phi(nodes,1);
            chi                               =  iso_coor(1);
            eta                               =  iso_coor(2);
            iota                              =  iso_coor(3);
            if isnumeric(str.f_e.symbolic_N_u)
                Nu                            =  (str.f_e.symbolic_N_u)';
            else
                Nu                            =  (eval(str.f_e.symbolic_N_u))';
            end
            if isnumeric(str.f_e.symbolic_N_F)
                NF                            =  (str.f_e.symbolic_N_F)';
            else
                NF                            =  (eval(str.f_e.symbolic_N_F))';
            end
            if isnumeric(str.f_e.symbolic_N_H)
                NH                            =  (str.f_e.symbolic_N_H)';
            else
                NH                            =  (eval(str.f_e.symbolic_N_H))';
            end
            if isnumeric(str.f_e.symbolic_N_J)
                NJ                            =  (str.f_e.symbolic_N_J)';
            else
                NJ                            =  (eval(str.f_e.symbolic_N_J))';
            end
            NSigmaF                           =  NF;
            NSigmaH                           =  NH;
            NSigmaJ                           =  NJ;
            DN_chi_u                          =  eval(str.f_e.symbolic_DN_chi_u);
            variables                         =  computing_variables(ielem_coarse,xelem,Xelem,phielem,Nu,DN_chi_u,NF,NH,NJ,NSigmaF,NSigmaH,NSigmaJ,str);
        end
    end
end