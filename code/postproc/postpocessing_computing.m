%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the postpocessing per element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  str                                                              =  postpocessing_computing(str)


switch str.data.actual_dim
    case {13,23}
          str                                                              =  postpocessing_computing_beams(str);
    otherwise
        
str.f_e                                                                    =  str.postproc.f_e;
switch str.data.shape
    case 1
str.f_e.N(:,[3 4 7 8])  =  str.f_e.N(:,[4 3 8 7]);
str.f_e.DN_chi(:,:,[3 4 7 8])  =  str.f_e.DN_chi(:,:,[4 3 8 7]);
end
 
str.quadrature.Chi                                                         =  str.nodal_positions_postprocessing_mesh;

%--------------------------------------------------------------------------
% 3. This function stores the different component of the stresses in vector
%    form.
%--------------------------------------------------------------------------
f                                                                          =  @(A,i,j) reshape(squeeze(A(i,j,:)),1,size(str.quadrature.Chi,1));
switch str.data.dim
    case 2
        g                                                                  =  @(B) [f(B,1,1);f(B,2,2);f(B,1,2)];
        dim                                                                =  2;
    case 3
        g                                                                  =  @(B) [f(B,1,1);f(B,2,2);f(B,3,3);f(B,1,2);f(B,1,3);f(B,2,3)];
        dim                                                                =  3;
end

%-----------------------------------------------------------
% 4. Get the variables to plot in the nodes. 
%-----------------------------------------------------------
%str.nodes_counter                                                         =  zeros(1,str.postproc.n_nodes);
for ielem=1:str.n_elem
    str.ielem                                                              =  ielem;
    nodes_elem_plot                                                        =  str.postproc.connectivity(ielem,:);
     switch str.data.shape
        case 1
            switch str.data.dim 
                case 2
                     nodes_elem_plot(:,[3 4])                              =  nodes_elem_plot(:,[4 3]); 
                case 3
                     nodes_elem_plot(:,[3 4 7 8])                          =  nodes_elem_plot(:,[4 3 8 7]); 
            end
         case 0
                      %nodes_elem_plot(:,[1 2 3 4])                         =  nodes_elem_plot(:,[2  4  1  3]);
                      %nodes_elem_plot(:,[1 2 3 4])                         =  nodes_elem_plot(:,[1  4  3  2]);
     end
%     end
    %nodes_elem                                                             =  str.postproc.connectivity(ielem,:);
    solution_nodes_elem                                                    =  str.connectivity(ielem,:);
    %solution_nodes_elem                                                    =  str.postproc.connectivity(ielem,:);
    str.properties.material_identifier                                     =  str.properties.material(ielem,1);
    %-----------------------------
    % Current time step.
    %-----------------------------
    gauss_level_information                                                =  postprocessing_element_nodes(str);    
    Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem);        
    xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem);
    phielem                                                                =  str.phi(solution_nodes_elem,1);        
    
    %Xelem_postproc                                                         =  gauss_level_information.Lagrangian_X;
    xelem_postproc                                                         =  gauss_level_information.Eulerian_x;
    phielem_postproc                                                       =  gauss_level_information.phi;    
    velocity_postproc                                                      =  gauss_level_information.velocity;
    acceleration_postproc                                                  =  gauss_level_information.acceleration;
    switch str.data.formulation.mixed_type
        case {'twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v3_mixed_formulation'}
             str                                                           =  gradients_micropolar_J(xelem,Xelem,str);
        case {'twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
             str                                                           =  gradients_micropolar_J_v4(xelem,Xelem,zeros(2,size(str.f_e.DN_chi_J,2)),str);
        otherwise
             str                                                           =  gradients(xelem,Xelem,phielem,str);
    end 
    %str                                                                    =  gradients(xelem_postproc,Xelem_postproc,phielem_postproc,str);
    %-----------------------------
    % x and phi postprocessing nodes. 
    %-----------------------------        
    str.postproc.Eulerian_x(:,nodes_elem_plot)                             =  xelem_postproc;
    str.postproc.phi(nodes_elem_plot)                                      =  phielem_postproc';
    str.postproc.velocity(:,nodes_elem_plot)                               =  velocity_postproc;
    str.postproc.acceleration(:,nodes_elem_plot)                           =  acceleration_postproc;
    %----------------------------------------
    % Initialising everything to zero
    %----------------------------------------
    F                                                                      =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    H                                                                      =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    J                                                                      =  zeros(size(str.postproc.f_e.N,2),1);
    SigmaF                                                                 =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    SigmaH                                                                 =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    SigmaJ                                                                 =  zeros(size(str.postproc.f_e.N,2),1);
    E0                                                                     =  zeros(dim,size(str.postproc.f_e.N,2));
    E                                                                      =  zeros(dim,size(str.postproc.f_e.N,2));
    D0                                                                     =  zeros(dim,size(str.postproc.f_e.N,2));
    S                                                                      =  zeros(dim,dim,size(str.postproc.f_e.N,2));    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Computation of variables for the vacuum and solid meshes (Only in BEM)
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    if str.vacuum.flag
       str.f_e.N                                                           =  str.postproc.f_e.N;
       str.f_e.DN_chi                                                      =  str.postproc.f_e.DN_chi;
       new_str                                                             =  gradients(xelem,Xelem,phielem,str);
       Jn                                                                  =  new_str.grad.J;
       Fn                                                                  =  new_str.grad.F;
       J                                                                   =  new_str.grad.J;
       F                                                                   =  new_str.grad.F;
       H                                                                   =  new_str.grad.H;
       E0                                                                  =  new_str.grad.E0;
       E                                                                   =  new_str.grad.E;
       Kirchhoff                                                           =  zeros(str.data.dim,str.data.dim,size(str.f_e.N,2));
        if ielem>str.vacuum_solid_mesh.postproc.last_solid_element
           D                                                               =  str.vacuum_permittivity*E;
           for igauss=1:size(str.f_e.N,2)
               D0(:,igauss)                                                =  str.vacuum_permittivity*inv(H(:,:,igauss))*D(:,igauss); 
           end            
        else            
           switch str.data.formulation_type          
                  %--------------------------------------------------------
                  %  Displacement potential formulation 
                  %--------------------------------------------------------
                  case 'displacement_potential_formulation' 
                       switch str.formulation
                             case 'D_to_E'
                                   [D0,str]                                =  internal_Newton_Raphson(str);
                             case 'E'
                                   D0                                      =  electric_displacement(str);
                       end
                  case 'mixed_formulation'                                 
                       gauss_level_information                             =  gauss_level_information_mixed_formulations(str);
                       switch str.data.formulation.mixed_type
                              case 'displacement_potential_mixed_formulation'
                                   D0                                      =  -gauss_level_information.DGammaDE0;
                              case 'full_mixed_formulation_electroelasticity_F_D0'
                                   D0                                      =  gauss_level_information.D0;                                                  
                              case 'full_mixed_formulation_electroelasticity_F_E0'
                                   D0                                      =  gauss_level_information.D0;                         
                              case 'full_mixed_formulation_electroelasticity_S_E0', 
                                   D0                                      =  gauss_level_information.D0;                         
                             case 'full_mixed_formulation_electroelasticity_S_D0'
                                   D0                                      =  gauss_level_information.D0;                         
                             case 'u_p_phi_compressible'
                                   D0                                      =  gauss_level_information.DGammaDE0;
                       end
           end
        end
    else
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Computation of variables when the vacuum is not present in the
    % calculations.
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------        
    option                                                                 =  0;
    switch str.data.formulation_type          
        %------------------------------------------------------------------
        %  Displacement potential formulation
        %------------------------------------------------------------------
        case 'displacement_potential_formulation' 
             oldstr                                                        =  str;
             str.f_e.N                                                     =  str.postproc.f_e.N;
             str.f_e.DN_chi                                                =  str.postproc.f_e.DN_chi;
             new_str                                                       =  gradients(xelem,Xelem,phielem,str);
             str                                                           =  oldstr;
             Jn                                                            =  new_str.grad.J;
             Fn                                                            =  new_str.grad.F;
             J                                                             =  new_str.grad.J;
             F                                                             =  new_str.grad.F;
             H                                                             =  new_str.grad.H;
             E0                                                            =  new_str.grad.E0;
             E                                                             =  new_str.grad.E;
             %-----------------------------------------
             % Second Piola and different component.
             %-----------------------------------------
             switch str.vacuum.flag
                 case 0
                     switch str.formulation
                         case 'D_to_E'
                             [S]                                           =  second_piola_d(str);
                         case 'E'
                             [S]                                           =  second_piola(str);
                     end
                 case 1
                     vacuum_identifier                                     =  max(str.properties.material);
                     switch str.properties.material(ielem)
                         case vacuum_identifier
                             str.formulation                               =  'E';
                             [S]                                           =  second_piola_vacuum(str);
                         otherwise
                             str.formulation                               =  str.old_formulation;
                             [S]                                           =  second_piola_d(str);
                     end
             end
             Kirchhoff                                                     =  S;
             %-----------------------------------------
             % Lagrangian electric displacement.
             %-----------------------------------------
             str.ielem                                                     =  ielem;
             switch str.vacuum.flag
                 case 0
                     switch str.formulation
                         case 'D_to_E'
                             [D0,str]                                      =  internal_Newton_Raphson(str);
                         case 'E'
                             D0                                            =  electric_displacement(str);
                     end
                 case 1
                     vacuum_identifier                                     =  max(str.properties.material);
                     switch str.properties.material(ielem)
                         case vacuum_identifier
                             str.formulation                               =  'E';
                             [D0]                                          =  electric_displacement_vacuum(str);
                         otherwise
                             str.formulation                               =  str.old_formulation;
                             [D0,str]                                      =  internal_Newton_Raphson(str);
                     end
             end
        %------------------------------------------------------------------
        %  Mixed formulations
        %------------------------------------------------------------------    
        case 'mixed_formulation'            
             str.properties.stabilisation_parameter_momentum=0;
             gauss_level_information                                       =  gauss_level_information_mixed_formulations(str);
             %gauss_level_information                                       =  gauss_level_information_mixed_formulations_final(str);
             switch str.data.formulation.mixed_type 
                    case 'displacement_formulation' 
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         F3D                                               =  F;
                         H3D                                               =  H;
                         switch str.data.dim
                             case 2
                                  F3D(3,3,:)                               =  ones(1,size(F3D,3));
                                  H3D(3,3,:)                               =  J;
                         end
                         for igauss=1:size(str.f_e.N,2)
                             DWDF                                          =  gauss_level_information.DWDF(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DWDF                                  =  [DWDF(1:2,1); 0;DWDF(3:4,1); zeros(3,1);DWDF(5,1)];
                             end
                             DWDF                                          =  reshape(DWDF,3,3);
                             DWDF                                          =  DWDF';
                             Piola1                                        =  DWDF;                             
                             DWDH                                          =  gauss_level_information.DWDH(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DWDH                                  =  [DWDH(1:2,1); 0;DWDH(3:4,1); zeros(3,1);DWDH(5,1)];
                             end
                             DWDH                                          =  reshape(DWDH,3,3);
                             DWDH                                          =  DWDH';
                             Piola2                                        =  Javier_double_cross_product(DWDH,F3D,1,1,3);                             
                             Piola3                                        =  gauss_level_information.DWDJ(igauss)*H3D(:,:,igauss);                             
                             Piola                                         =  Piola1 + Piola2 + Piola3;                         
                             Kirchhoff(:,:,igauss)                         =  Piola; 
                             SigmaF(:,:,igauss)                            =  Piola1;
                             SigmaH(:,:,igauss)                            =  DWDH;
                             SigmaJ(igauss)                                =  gauss_level_information.DWDJ(igauss);
                         end 
                    case 'full_mixed_formulation' 
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                    case 'F_J_mixed_formulation' 
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                    case 'u_J_mixed_formulation' 
                         Fn                                                =  str.grad.F;
                         Hn                                                =  str.grad.H;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.F;
                         J                                                 =  gauss_level_information.J';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                    case 'H_J_mixed_formulation' 
                         oldstr                                            =  str;
                         str.f_e.N                                         =  str.postproc.f_e.N;
                         str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
                         str                                               =  gradients(xelem,Xelem,phielem,str);
                         Fn                                                =  str.grad.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.DWDF;
                         SSSigmaF                                          =  zeros(3,3,size(str.quadrature.Chi,1));
                         for igauss=1:size(str.quadrature.Chi,1)
                             switch str.data.dim
                                 case 2
                                     SSigmaF                               =  [SigmaF(1:2,igauss); 0; SigmaF(3:4,igauss); zeros(3,1);SigmaF(5,igauss)];
                                 case 3 
                                     SSigmaF                               =  SigmaF(:,igauss);
                             end
                             SSigmaF                                       =  reshape(SSigmaF,3,3);
                             SSigmaF                                       =  SSigmaF';
                             SSSigmaF(:,:,igauss)                          =  SSigmaF;
                         end
                         SigmaF                                            =  SSSigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         str                                               =  oldstr;                         
                         Kirchhoff                                         =  Piola;
                    case {'complementary_energy','stabilised_complementary_energy','stabilised_complementary_energy_1','stabilised_complementary_energy_2','stabilised_complementary_energy_3'}
                         oldstr                                            =  str;
                         str.f_e.N                                         =  str.postproc.f_e.N;
                         str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
                         str                                               =  gradients(xelem,Xelem,phielem,str);
                         provF                                             =  gauss_level_information.DGDSF;
                         provH                                             =  gauss_level_information.DGDSH;
                         provJ                                             =  gauss_level_information.DGDSJ;
                         for igauss=1:size(str.quadrature.Chi,1)
                             Fn(:,:,igauss)                                =  (reshape(provF(:,igauss),3,3))';
                             Hn(:,:,igauss)                                =  (reshape(provH(:,igauss),3,3))';
                         end
                         Jn                                                =  provJ;
                         F                                                 =  Fn;
                         H                                                 =  Hn;
                         J                                                 =  Jn;
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         str                                               =  oldstr;                         
                         Kirchhoff                                         =  Piola; 
                    case {'full_incompressible','compressible_two_field'}
                         oldstr                                            =  str;
                         str.f_e.N                                         =  str.postproc.f_e.N;
                         str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
                         str                                               =  gradients(xelem,Xelem,phielem,str);
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         PPiola                                            =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
                         for igauss=1:size(str.quadrature.Chi,1)
                             Piola(:,:,igauss)                             =  PPiola(:,:,igauss) + 0*gauss_level_information.p(igauss)*eye(3);                         
                         end                                                  
                         Kirchhoff                                         =  Piola;
                         str                                               =  oldstr;
                         p_old                                             =  gauss_level_information.p;
                    case {'nearly_incompressible_three_field'}
                         oldstr                                            =  str;
                         str.f_e.N                                         =  str.postproc.f_e.N;
                         str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
                         str                                               =  gradients(xelem,Xelem,phielem,str);
                         Fn                                                =  str.grad.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         PPiola                                            =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
                         for igauss=1:size(str.quadrature.Chi,1)
                             Piola(:,:,igauss)                             =  PPiola(:,:,igauss) + 0*gauss_level_information.p(igauss)*eye(3);                         
                         end                                                  
                         Kirchhoff                                         =  Piola;
                         str                                               =  oldstr;
                         p_old                                             =  gauss_level_information.p;
                    case {'nearly_incompressible_two_field'}
                         oldstr                                            =  str;
                         str.f_e.N                                         =  str.postproc.f_e.N;
                         str.f_e.DN_chi                                    =  str.postproc.f_e.DN_chi;
                         str                                               =  gradients(xelem,Xelem,phielem,str);
                         Fn                                                =  str.grad.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         PPiola                                            =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);
                         for igauss=1:size(str.quadrature.Chi,1)
                             Piola(:,:,igauss)                             =  PPiola(:,:,igauss) + 0*gauss_level_information.p(igauss)*eye(3);                         
                         end                                                  
                         Kirchhoff                                         =  Piola;
                         str                                               =  oldstr;
                         p_old                                             =  gauss_level_information.p;
                    case {'superfull_mixed_incompressible_formulation','full_mixed_incompressible_formulation'}
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         p_old                                             =  gauss_level_information.p;
                    case 'displacement_potential_mixed_formulation'
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         E0                                                =  str.grad.E0;
                         F3D                                               =  F;
                         H3D                                               =  H;
                         switch str.data.dim
                             case 2
                                  F3D(3,3,:)                               =  ones(1,size(F3D,3));
                                  H3D(3,3,:)                               =  J;
                         end
                         for igauss=1:size(str.f_e.N,2)
                             DGammaDF                                      =  gauss_level_information.DGammaDF(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DGammaDF                              =  [DGammaDF(1:2,1); 0;DGammaDF(3:4,1); zeros(3,1);DGammaDF(5,1)];
                             end
                             DGammaDF                                      =  reshape(DGammaDF,3,3);
                             DGammaDF                                      =  DGammaDF';
                             Piola1                                        =  DGammaDF;                             
                             DGammaDH                                      =  gauss_level_information.DGammaDH(:,igauss);
                             switch str.data.dim  
                                 case 2 
                                     DGammaDH                              =  [DGammaDH(1:2,1); 0;DGammaDH(3:4,1); zeros(3,1);DGammaDH(5,1)];
                             end
                             DGammaDH                                      =  reshape(DGammaDH,3,3);
                             DGammaDH                                      =  DGammaDH';
                             Piola2                                        =  Javier_double_cross_product(DGammaDH,F3D,1,1,3);                             
                             Piola3                                        =  gauss_level_information.DGammaDJ(igauss)*H3D(:,:,igauss);                             
                             Piola                                         =  Piola1 + Piola2 + Piola3;                         
                             Kirchhoff(:,:,igauss)                         =  Piola; 
                             SigmaF(:,:,igauss)                            =  Piola1;
                             SigmaH(:,:,igauss)                            =  DGammaDH;
                             SigmaJ(igauss)                                =  gauss_level_information.DGammaDJ(igauss);
                             D0(:,igauss)                                  =  -gauss_level_information.DGammaDE0(:,igauss);
                         end  
                    case 'u_p_phi_compressible'  
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         E0                                                =  str.grad.E0;
                         F3D                                               =  F;
                         H3D                                               =  H;
                         switch str.data.dim
                             case 2
                                  F3D(3,3,:)                               =  ones(1,size(F3D,3));
                                  H3D(3,3,:)                               =  J;
                         end
                         for igauss=1:size(str.f_e.N,2)
                             DGammaDF                                      =  gauss_level_information.DGammaDF(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DGammaDF                              =  [DGammaDF(1:2,1); 0;DGammaDF(3:4,1); zeros(3,1);DGammaDF(5,1)];
                             end
                             DGammaDF                                      =  reshape(DGammaDF,3,3);
                             DGammaDF                                      =  DGammaDF';
                             Piola1                                        =  DGammaDF;                             
                             DGammaDH                                      =  gauss_level_information.DGammaDH(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DGammaDH                              =  [DGammaDH(1:2,1); 0;DGammaDH(3:4,1); zeros(3,1);DGammaDH(5,1)];
                             end
                             DGammaDH                                      =  reshape(DGammaDH,3,3);
                             DGammaDH                                      =  DGammaDH';
                             Piola2                                        =  Javier_double_cross_product(DGammaDH,F3D,1,1,3);                             
                             Piola3                                        =  gauss_level_information.DGammaDJ(igauss)*H3D(:,:,igauss);                             
                             Piola                                         =  Piola1 + Piola2 + Piola3;                         
                             Kirchhoff(:,:,igauss)                         =  Piola; 
                             SigmaF(:,:,igauss)                            =  Piola1;
                             SigmaH(:,:,igauss)                            =  DGammaDH;
                             SigmaJ(igauss)                                =  gauss_level_information.DGammaDJ(igauss);
                             D0(:,igauss)                                  =  -gauss_level_information.DGammaDE0(:,igauss);
                             p_old                                         =  gauss_level_information.p;
                         end 
                    case 'u_p_phi_incompressible'  
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         E0                                                =  str.grad.E0;
                         F3D                                               =  F;
                         H3D                                               =  H;
                         switch str.data.dim
                             case 2
                                  F3D(3,3,:)                               =  ones(1,size(F3D,3));
                                  H3D(3,3,:)                               =  J;
                         end
                         for igauss=1:size(str.f_e.N,2)
                             DGammaDF                                      =  gauss_level_information.DGammaDF(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DGammaDF                              =  [DGammaDF(1:2,1); 0;DGammaDF(3:4,1); zeros(3,1);DGammaDF(5,1)];
                             end
                             DGammaDF                                      =  reshape(DGammaDF,3,3);
                             DGammaDF                                      =  DGammaDF';
                             Piola1                                        =  DGammaDF;                             
                             DGammaDH                                      =  gauss_level_information.DGammaDH(:,igauss);
                             switch str.data.dim
                                 case 2
                                     DGammaDH                              =  [DGammaDH(1:2,1); 0;DGammaDH(3:4,1); zeros(3,1);DGammaDH(5,1)];
                             end
                             DGammaDH                                      =  reshape(DGammaDH,3,3);
                             DGammaDH                                      =  DGammaDH';
                             Piola2                                        =  Javier_double_cross_product(DGammaDH,F3D,1,1,3);                             
                             Piola3                                        =  gauss_level_information.DGammaDJ(igauss)*H3D(:,:,igauss);                             
                             Piola                                         =  Piola1 + Piola2 + Piola3;                         
                             Kirchhoff(:,:,igauss)                         =  Piola + gauss_level_information.P(igauss)*H3D(:,:,igauss);                             
                             SigmaF(:,:,igauss)                            =  Piola1;
                             SigmaH(:,:,igauss)                            =  DGammaDH;
                             SigmaJ(igauss)                                =  gauss_level_information.DGammaDJ(igauss);
                             D0(:,igauss)                                  =  -gauss_level_information.DGammaDE0(:,igauss);
                             p_old                                         =  gauss_level_information.P;
                         end 
                    case 'full_mixed_formulation_electroelasticity_F_D0'
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         E0                                                =  gauss_level_information.DUDD0;
                         D0                                                =  gauss_level_information.D0;                                                  
                    case 'full_mixed_formulation_electroelasticity_F_E0'
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         E0                                                =  gauss_level_information.E0;
                         D0                                                =  gauss_level_information.D0;                         
                    case 'full_mixed_formulation_electroelasticity_S_E0', 
                         Fn                                                =  reshape(gauss_level_information.DUpsilonDSF',str.data.dim,str.data.dim,size(gauss_level_information.DUpsilonDSF,1)*size(gauss_level_information.DUpsilonDSF,2)/str.data.dim^2);
                         Jn                                                =  gauss_level_information.DUpsilonDSJ;
                         F                                                 =  reshape(gauss_level_information.DUpsilonDSF',str.data.dim,str.data.dim,size(gauss_level_information.DUpsilonDSF,1)*size(gauss_level_information.DUpsilonDSF,2)/str.data.dim^2);
                         H                                                 =  reshape(gauss_level_information.DUpsilonDSH',str.data.dim,str.data.dim,size(gauss_level_information.DUpsilonDSH,1)*size(gauss_level_information.DUpsilonDSH,2)/str.data.dim^2);
                         J                                                 =  gauss_level_information.DUpsilonDSJ;
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         E0                                                =  gauss_level_information.E0;
                         D0                                                =  gauss_level_information.D0;                         
                    case 'full_mixed_formulation_electroelasticity_S_D0'
                         Jn                                                =  gauss_level_information.DSigmaDSJ;
                         J                                                 =  gauss_level_information.DSigmaDSJ;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         for igauss=1:size(str.quadrature.Chi,1)                        
                         Fn(:,:,igauss)                                    =  (reshape(gauss_level_information.DSigmaDSF(:,igauss),str.data.dim,str.data.dim,1))';
                         F(:,:,igauss)                                     =  (reshape(gauss_level_information.DSigmaDSF(:,igauss),str.data.dim,str.data.dim,1))';
                         H(:,:,igauss)                                     =  (reshape(gauss_level_information.DSigmaDSH(:,igauss),str.data.dim,str.data.dim,1))';
                         end 
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         E0                                                =  -gauss_level_information.DSigmaDD0;
                         D0                                                =  gauss_level_information.D0;                         
                         option                                            =  SigmaV;
                    case 'full_mixed_formulation_electroelasticity_F_D0_V'  
                          Fx =  str.grad.F; 
                          Hx =  str.grad.H;
                          Jx =  str.grad.J;
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                          %F                                                 =  str.grad.F;
                          %H                                                 =  str.grad.H;
                          %J                                                 =  str.grad.J;
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         SigmaD0                                           =  gauss_level_information.DUDD0;
                         SigmaV                                            =  gauss_level_information.DUDV;
                         D0                                                =  gauss_level_information.D0;                                                  
                         V                                                 =  gauss_level_information.V;                                                  
                         for igauss=1:size(str.quadrature.Chi,1)                        
                         E0(:,igauss)                                      =  SigmaD0(:,igauss) + Fn(:,:,igauss)'*SigmaV(:,igauss);
                         end
                         option  =  cell(1,3); 
                         option{1,1}                                       =  SigmaV;
                         option{1,2}                                       =  D0;
                         option{1,3}                                       =  V;
                    case 'full_mixed_formulation_electroelasticity_F_E0_V'
                         Fn                                                =  gauss_level_information.F;
                         Jn                                                =  gauss_level_information.J;
                         F                                                 =  gauss_level_information.F;
                         H                                                 =  gauss_level_information.H;
                         J                                                 =  gauss_level_information.J';
                         SigmaF                                            =  gauss_level_information.SigmaF;
                         SigmaH                                            =  gauss_level_information.SigmaH;
                         SigmaD0                                           =  gauss_level_information.SigmaD0;
                         SigmaV                                            =  gauss_level_information.SigmaV;
                         SigmaJ                                            =  gauss_level_information.SigmaJ';
                         Piola                                             =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);                         
                         Kirchhoff                                         =  Piola;
                         D0                                                =  gauss_level_information.D0;                         
                         for igauss=1:size(str.quadrature.Chi,1)                        
                         E0(:,igauss)                                      =  SigmaD0(igauss) + Fn(:,:,igauss)*SigmaV(:,igauss);
                         end
                         option  =  cell(1,3);
                         option{1,1}                                       =  SigmaV;
                         option{1,2}                                       =  D0;
                    case {'twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v3_mixed_formulation','twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
                         Fn                                                =  str.grad.F;
                         Jn                                                =  str.grad.J;
                         F                                                 =  str.grad.F;
                         H                                                 =  str.grad.H;
                         J                                                 =  str.grad.J;
                         gradJ                                             =  str.grad.DJ;
                         ratioJ                                            =  min(str.grad.J)/max(str.grad.J);
                         ratioJ                                            =  ratioJ*ones(size(J,1),size(J,2));
                         for igauss=1:size(str.f_e.N,2)
                             DWDF                                          =  gauss_level_information.DWDF(:,:,igauss);
                             Piola1                                        =  DWDF;                             
                             Piola3                                        =  gauss_level_information.DWDJ(igauss)*H(:,:,igauss);                             
                             Piola                                         =  Piola1 + Piola3;                         
                             Kirchhoff(:,:,igauss)                         =  Piola; 
                             SigmaF(:,:,igauss)                            =  Piola1;
                             SigmaH(:,:,igauss)                            =  zeros(2,2);
                             SigmaJ(igauss)                                =  gauss_level_information.DWDJ(igauss);
                         end 
             end 
    
    end
    end
    %----------------------------------------------------------------------
    % Compute norms of tensors
    %----------------------------------------------------------------------    
    FL2norm                                                                =  zeros(size(str.f_e.N,2),1);
    for igauss=1:size(str.f_e.N,2)
        FL2norm(igauss)                                                    =  sqrt(trace(F(:,:,igauss)'*F(:,:,igauss)));
    end
    
    %----------------------------------------------------------------------
    % Eigenvalues of the second derivative of the energy functional with
    % respect to the deformation gradient F. This is done for stability
    % analysis purposes.
    %---------------------------------------------------------------------- 
    eigenvalues                                                            =  zeros(str.data.dim^2,size(str.quadrature.Chi,1));
    if 0
    for igauss=1:size(str.quadrature.Chi,1)
        F_                                                                 =  F(:,:,igauss);
        H_                                                                 =  H(:,:,igauss);
        J_                                                                 =  J(igauss);
        eigenvalues(:,igauss)                                              =  eigenvalues_constitutive_tensor(str,F_,H_,J_);
    end    
    end
    %-----------------------------------------
    % Cauchy stress and different component.
    %-----------------------------------------
    sigma                                                                  =  cauchy_stress(str,Kirchhoff,F,H,J,SigmaF,SigmaH,SigmaJ,option,gauss_level_information);
    %-----------------------------------------
    % Volumetric and deviatoric decomposition
    % of Cauchy and Piola.
    %-----------------------------------------
    [sigma_deviatoric,sigma_volumetric,...
        S_deviatoric,S_volumetric,p]                                       =  volumetric_deviatoric_cauchy_Piola_decomposition(str,sigma);    
    switch str.data.formulation.mixed_type
        %-------------------------------------------------------------
        % Update pressure
        %-------------------------------------------------------------
        case {'full_incompressible','compressible_two_field','nearly_incompressible_three_field','nearly_incompressible_two_field',...
                'full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation','u_p_phi_compressible'}
            
             p                                                             =  p_old';
    end    
    %-----------------------------------------            
    % Eulerian electric displacement and different components.
    %-----------------------------------------
    [eul_D,E,E_norm]                                                       =  elec_disp(str,D0,E0);
    %-----------------------------------------            
    % Eulerian electric displacement and different components.
    %-----------------------------------------
    factor                                                                 =  electric_breakdown_factor(str,E_norm,F);
    %----------------------------------------------------------------
    % Second Piola and different contributions. 
    %----------------------------------------------------------------
    str.postproc.stress(:,nodes_elem_plot)                                 =  str.postproc.stress(:,nodes_elem_plot) + g(S);
    str.postproc.stress_deviatoric(:,nodes_elem_plot)                      =  str.postproc.stress_deviatoric(:,nodes_elem_plot) + g(S_deviatoric);
    str.postproc.stress_volumetric(:,nodes_elem_plot)                      =  str.postproc.stress_volumetric(:,nodes_elem_plot) + g(S_volumetric);
    %----------------------------------------------------------------
    % First Piola stress conjugates.
    %----------------------------------------------------------------
    str.postproc.SigmaF(:,:,nodes_elem_plot)                               =  str.postproc.SigmaF(:,:,nodes_elem_plot) + twoD_treatment(SigmaF);
    str.postproc.SigmaH(:,:,nodes_elem_plot)                               =  str.postproc.SigmaH(:,:,nodes_elem_plot) + twoD_treatment(SigmaH);
    str.postproc.SigmaJ(nodes_elem_plot,1)                                 =  str.postproc.SigmaJ(nodes_elem_plot,1) + SigmaJ;
    str.postproc.First_Piola(:,:,nodes_elem_plot)                          =  str.postproc.First_Piola(:,:,nodes_elem_plot) + twoD_treatment(Kirchhoff);
    %----------------------------------------------------------------
    % Cauchy stress and different contributions.
    %----------------------------------------------------------------
    str.postproc.sigma(:,nodes_elem_plot)                                  =  str.postproc.sigma(:,nodes_elem_plot) +  g(sigma);
    str.postproc.sigma_deviatoric(:,nodes_elem_plot)                       =  str.postproc.sigma_deviatoric(:,nodes_elem_plot) + g(sigma_deviatoric);
    str.postproc.sigma_volumetric(:,nodes_elem_plot)                       =  str.postproc.sigma_volumetric(:,nodes_elem_plot) + g(sigma_volumetric);
    str.postproc.sigma_pressure(nodes_elem_plot,1)                         =  str.postproc.sigma_pressure(nodes_elem_plot,1) + p;
    %----------------------------------------------------------------
    % Electric displacement. Eulerian part.
    %----------------------------------------------------------------
    str.postproc.D(:,nodes_elem_plot)                                      =  str.postproc.D(:,nodes_elem_plot) + eul_D;
    %----------------------------------------------------------------
    % Electric displacement. Lagrangian part.
    %----------------------------------------------------------------
    str.postproc.D0(:,nodes_elem_plot)                                     =  str.postproc.D0(:,nodes_elem_plot) + D0;
    %----------------------------------------------------------------
    % Electric field.
    %----------------------------------------------------------------
    str.postproc.E0(:,nodes_elem_plot)                                     =  str.postproc.E0(:,nodes_elem_plot) + E0;
    str.postproc.E(:,nodes_elem_plot)                                      =  str.postproc.E(:,nodes_elem_plot) + E;
    str.postproc.E_norm(nodes_elem_plot,1)                                 =  str.postproc.E_norm(nodes_elem_plot,1) + E_norm;
    str.postproc.electric_breakdown_factor(nodes_elem_plot,1)              =  str.postproc.electric_breakdown_factor(nodes_elem_plot,1) + factor;
    %----------------------------------------------------------------
    % Deformation.
    %----------------------------------------------------------------
    str.postproc.F(:,:,nodes_elem_plot)                                    =  str.postproc.F(:,:,nodes_elem_plot) + twoD_treatment(F);
    str.postproc.FL2norm(nodes_elem_plot)                                  =  str.postproc.FL2norm(nodes_elem_plot) + FL2norm;                  
    str.postproc.H(:,:,nodes_elem_plot)                                    =  str.postproc.H(:,:,nodes_elem_plot) + twoD_treatment(H);
    str.postproc.J(nodes_elem_plot,1)                                      =  str.postproc.J(nodes_elem_plot,1) + J;
    switch str.data.formulation.mixed_type
        case {'twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v3_mixed_formulation','twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
             str.postproc.ratioJ(nodes_elem_plot,1)                        =  str.postproc.ratioJ(nodes_elem_plot,1) + ratioJ;
             str.postproc.gradJ(:,nodes_elem_plot)                         =  str.postproc.gradJ(:,nodes_elem_plot) + gradJ;
                  
    end
    %----------------------------------------------------------------
    % Second derivative of the energy function with respect to F. 
    %----------------------------------------------------------------
    str.postproc.eigenvalues(:,nodes_elem_plot)                            =  str.postproc.eigenvalues(:,nodes_elem_plot) + eigenvalues;
    %----------------------------- 
    % Previous time step in order to get J.
    %-----------------------------
    if 0
    switch str.data.analysis
        case 'dynamic'
            if str.time_iteration>=2
                switch str.formulation
                    case 'E'
                        xelemn1                                            =  str.Eulerian_xn1(:,nodes_elem);
                        phielemn1                                          =  str.phin1(nodes_elem);
                        [str]                                              =  gradients(xelemn1,Xelem,phielemn1,str);
                        [D0n1]                                             =  electric_displacement(str);
                        [J_electric,J0_electric]                           =  Intensity_vector(str,D0,D0n1,Jn,Fn);
                        str.postproc.J_electric(:,...
                            nodes_elem_plot)                               =  str.postproc.J_electric(:,nodes_elem_plot) + J_electric;
                        str.postproc.J0_electric(:,...
                            nodes_elem_plot)                               =  str.postproc.J0_electric(:,nodes_elem_plot) + J0_electric;
                    case 'D_to_E'
                        xelemn1                                            =  str.Eulerian_xn1(:,nodes_elem);
                        phielemn1                                          =  str.phin1(nodes_elem);
                        [str]                                              =  gradients(xelemn1,Xelem,phielemn1,str);
                        [D0n1,str]                                         =  internal_Newton_Raphson(str);
                        [J_electric,J0_electric]                           =  Intensity_vector(str,D0,D0n1,Jn,Fn);
                        str.postproc.J_electric(:,...
                            nodes_elem_plot)                               =  str.postproc.J_electric(:,nodes_elem_plot) + J_electric;
                        str.postproc.J0_electric(:,...
                            nodes_elem_plot)                               =  str.postproc.J0_electric(:,nodes_elem_plot) + J0_electric;
                end
            end
    end
    end
    switch str.data.elastoplasticity
        case 1
            [S,Dv,new_str]                                                 =  Radial_Return_Mapping_algorithm(str);
            for inode=1:size(str.postproc.connectivity,2)
                lambda_p                                                   =  eig(str.solid.kinematic.plastic.C_p(:,:,ielem,inode));
                lambda_p                                                   =  sqrt(lambda_p);
                node                                                       =  str.postproc.connectivity(ielem,inode);
                str.postproc.hardening_variable(node)                      =  str.postproc.hardening_variable(inode) + new_str.solid.kinematic.plastic.hardening_variable(ielem,inode);
                str.postproc.lambda_p(:,node)                              =  str.postproc.lambda_p(:,node)  +  lambda_p;
            end
        case 0
            for inode=1:size(str.postproc.connectivity,2)
                lambda_p                                                   =  sqrt(eig(eye(dim)));
                node                                                       =  str.postproc.connectivity(ielem,inode);
                str.postproc.lambda_p(:,node)                              =  lambda_p;
                str.postproc.hardening_variable(node)                      =  0;
            end
    end    
    %----------------------------------------------------------------
    % Element contributions on each node.
    %----------------------------------------------------------------
    %str.nodes_counter(nodes_elem_plot)                                    =  str.nodes_counter(nodes_elem_plot) + 1;
end
    str.postproc.phi=0*ones(size(str.postproc.phi,1),1);

end
function tensor   =  twoD_treatment(tensor)
     switch str.data.dim
         case 2
             switch size(tensor,1)
                 case 3
                      tensor(3,:,:)                                        =  [];
                      tensor(:,3,:)                                        =  [];
             end
     end      
end

end