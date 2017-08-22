function [I,W_dot_elec]                 =  Electric_Intensity(str)

data                                    =  cell(3);

for iedge =1:size(str.surface_elem,2)
    element_info                        =  str.surface_elem{iedge};
    element                             =  element_info(1,1);
    size_surf                           =  size(element_info,2);
    str.aux.nodes                       =  element_info(2:size_surf-2);
    str.aux.global_nodes                =  str.connectivity(element,str.aux.nodes);
    %-------------------------------------------
    % Material associated.
    %-------------------------------------------
    str.properties.material_identifier  =  str.properties.material(element,1);
    %-------------------------------------------
    % Store old_str.
    %-------------------------------------------
    old_str                             =  str;
    %-------------------------------------------
    % Coordinates and electric potential for
    % current and previous time step.
    %-------------------------------------------
    xelem                               =  str.Eulerian_x(:,str.aux.global_nodes);
    Xelem                               =  str.Lagrangian_X(:,str.aux.global_nodes);
    phielem                             =  str.phi(str.aux.global_nodes);
    xelemn1                             =  str.Eulerian_xn1(:,str.aux.global_nodes);
    phielemn1                           =  str.phin1(str.aux.global_nodes);
    %------------------------------------------------------------------
    % I get the normal vector to the edge for the current time step.
    %------------------------------------------------------------------
    str.quadrature.Chi                  =  str.flux_calc.quadrature.Chi;
    str.quadrature.W_v                  =  str.flux_calc.quadrature.W_v;
    str.f_e.N                           =  str.flux_calc.f_e.N;
    str.f_e.DN_chi                      =  str.flux_calc.f_e.DN_chi;
    str.data.old_dim                    =  str.data.dim;
    str.data.dim                        =  100;
    str.n_node_elem                     =  size(str.aux.global_nodes,2);
    [str]                               =  gradients(xelem,Xelem,phielem,str);
    Normal_vector                       =  str.grad.Normal_iota;
    str                                 =  old_str;
    %------------------------------------------------------------------
    % Identification of the gauss points in the edge or surface element
    % with its position in the isoparametric system of the complete element.
    %------------------------------------------------------------------
    data{1}                             =  xelem;
    data{2}                             =  xelem;
    data{3}                             =  phielem;
    [Chi]                               =  mapping_igauss_surface_vol_elements(data,old_str,str);
    str.quadrature.Chi                  =  Chi;
    %--------------------------------------------------------------------------
    % D0 in current time step.
    %--------------------------------------------------------------------------
    [str]                               =  gradients(xelem,Xelem,phielem,str);
    [D0]                                =  electric_displacement(str);
    str                                 =  old_str;
    %--------------------------------------------------------------------------
    % D0 in previous time step.
    %--------------------------------------------------------------------------
    [str]                               =  gradients(xelemn1,Xelem,phielemn1,str);
    [D0n1]                              =  electric_displacement(str);
    str                                 =  old_str;
    %------------------------------------------------------------------
    % Intensity, energy power.
    %------------------------------------------------------------------
        for iloop1=1:length(str.quadrature.Chi)
            DX_chi                      =  str.grad.DX_chi(:,:,iloop1);
            J_t                         =  DX_chi';
            J_t                         =  abs(det(J_t));
            W                           =  str.quadrature.W_v(iloop1);
            %---------------------------------------------------------------------
            % Obtain Lagrangian electric displacement.
            %---------------------------------------------------------------------
            D0_                         =  D0(:,iloop1);
            D0n1_                       =  D0n1(:,iloop1);
            D0Dt                        =  (D0_ - D0n1_)/Dt;
            phi                         =  phielem*str.flux_calc.f_e.N(:,iloop1);  % Electric potential in the gauss point considered.
            %-------------------------------------------------------------------
            % Force vector (motion part).
            %-------------------------------------------------------------------
            I                           =  I  + D0Dt'*N*W*J_t;
            W_dot_elec                  =  W_dot_elec  + D0Dt'*Normal_vector*phi*W*J_t;
        end
end
    
end



function [Chi]                      =  mapping_igauss_surface_vol_elements(data,old_str,str)

xelem                               =  data{1};
Xelem                               =  data{2};
phielem                             =  data{3};
flux_dim                            =  str.data.dim - 1;
%------------------------------------------------------------------
% I identify the gauss points in the edge . 
%------------------------------------------------------------------
isopar_nodes                        =  str.nodal_positions(str.aux.nodes,:);
dot_product                         =  zeros(str.data.dim,1);
vector                              =  zeros(size(str.aux.global_nodes,2),str.data.dim);
for idim =1:str.data.dim
    for inode=1:size(str.aux.global_nodes,2)
        vector(inode,idim)          =  isopar_nodes(inode,idim);
    end
    dot_product(idim)               =  vector(:,idim)'*ones(size(str.aux.global_nodes,2),1);
    if abs(dot_product(idim)-size(str.aux.global_nodes,2))<1e-7
       coor_value                   =  vector(inode,idim);
       fixed_flux_dimension         =  [idim coor_value];
    end
end
%------------------------------------------------------------------
% Gauss points in the midle of the edge for the surface element. 
%------------------------------------------------------------------
str.quadrature.Chi                  =  zeros(1,str.data.dim);  % it has to be zero in the dimension(s) out of the element
str.quadrature.Chi(:,...
    fixed_flux_dimension(1))        =  fixed_flux_dimension(2); 
%--------------------------------------------------------------------------
% Tangent vectors in the middle of the edge (or face) for the complete 
% element (with all its nodes, not just with the nodes of the face).
%--------------------------------------------------------------------------
[str]                               =  gradients(xelem,Xelem,phielem,str);
N1                                  =  str.grad.DX_chi(:,1,:);
N2                                  =  str.grad.DX_chi(:,2,:);
N3                                  =  str.grad.DX_chi(:,3,:);
str                                 =  old_str;
%--------------------------------------------------------------------------
% Tangent vectors in the middle of the edge (or face) for the surface 
% element (just with the nodes of the face).
%--------------------------------------------------------------------------
str.quadrature.Chi                  =  zeros(1,flux_dim);
str.data.dim                        =  100;
str.f_e.N                           =  str.flux_calc.f_e.N;
str.f_e.DN_chi                      =  str.flux_calc.f_e.DN_chi;
str.n_node_elem                     =  size(str.aux.global_nodes,2);
[str]                               =  gradients(xelem,Xelem,phielem,str);
N1_edge                             =  str.grad.Normal_chi;
N2_edge                             =  str.grad.Normal_eta;
N3_edge                             =  str.grad.Normal_iota;
str                                 =  old_str;
%--------------------------------------------------------------------------
% Rotation matrix between the tangent vectors from the complete and surface
% (or line) elements. 
%--------------------------------------------------------------------------
Lambda_rot                          =  N1_edge*N1' + N2_edge*N2' + N3_edge*N3';
%--------------------------------------------------------------------------
% We obtain the gauss point of the surface element in the isoparametric
% coordinate system of the complete element.
%--------------------------------------------------------------------------
Chi                                 =  zeros(size(str.flux_calc.quadrature.Chi,1),str.data.dim);
Chi(:,1:flux_dim)                   =  str.flux_calc.quadrature.Chi;
Chi                                 =  Lambda_rot*Chi;
end

end












