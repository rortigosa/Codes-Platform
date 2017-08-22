function [str]                                                 =  former_postprocessing(str)


%--------------------------------------------------------------------------
% 1. Finite element technology and Gauss integration for postprocessing
%--------------------------------------------------------------------------
% if ~isfield(str,'postproc')
 str.postproc.output_type                                      =  'averaging';
% end
%--------------------------------------------------------------------------
% 1.1 Remeshing of the domain in
% order to have just linear elements. 
%--------------------------------------------------------------------------
%[str]                                                         =  linear_multielement_remeshing(str) 
%--------------------------------------------------------------------------
% 1.2 Reordering of shape functions and connectivities.
%--------------------------------------------------------------------------
switch str.data.dim
    case {13,23}
        str.Eulerian_x                                         =  str.solid.BEAM_SHELL.continuum.Eulerian_x;
end
switch str.data.dim
    case {2,3}
         old_str.data.dim                                      =  str.data.dim;
    case {13,23}
         old_str.data.dim                                      =  str.data.dim;
         str.data.dim                                          =  3;         
end
%--------------------------------------------------------------------------
% Evaluation of shape functions in the nodes of the element with linear
% degree of approximation.
%--------------------------------------------------------------------------
str                                                            =  isoparametric_nodes(str);
str                                                            =  shape_function_displacement(str);
switch str.data.formulation_type
    case 'mixed_formulation'
         str                                                   =  shape_function_nodes_mixed_formulation(str);
end

N                                                              =  str.f_e.N;
DN                                                             =  str.f_e.DN_chi;
connectivity                                                   =  str.connectivity;

%--------------------------------------------------------------------------
%  Reordering of the shape functions in agreement with our numbering
%  criterion.
%--------------------------------------------------------------------------
switch str.data.shape
    case 1
        switch str.data.dim
            case 2
                N(3,:)                                         =  str.f_e.N(4,:);
                N(4,:)                                         =  str.f_e.N(3,:);
                DN(:,3,:)                                      =  str.f_e.DN_chi(:,4,:);
                DN(:,4,:)                                      =  str.f_e.DN_chi(:,3,:);
                connectivity(:,3)                              =  str.connectivity(:,4);
                connectivity(:,4)                              =  str.connectivity(:,3);
            case 3
                N(3,:)                                         =  str.f_e.N(4,:);
                N(4,:)                                         =  str.f_e.N(3,:);
                N(7,:)                                         =  str.f_e.N(8,:);
                N(8,:)                                         =  str.f_e.N(7,:);
                DN(:,3,:)                                      =  str.f_e.DN_chi(:,4,:);
                DN(:,4,:)                                      =  str.f_e.DN_chi(:,3,:);
                DN(:,7,:)                                      =  str.f_e.DN_chi(:,8,:);
                DN(:,8,:)                                      =  str.f_e.DN_chi(:,7,:);
                connectivity(:,3)                              =  str.connectivity(:,4);
                connectivity(:,4)                              =  str.connectivity(:,3);
                connectivity(:,7)                              =  str.connectivity(:,8);
                connectivity(:,8)                              =  str.connectivity(:,7);
        end
end
str.f_e.N                                                      =  N;
str.f_e.DN_chi                                                 =  DN;
str.connectivity                                               =  connectivity;
         
         
switch str.data.example_comparison
    case 0
         %--------------------------------------------------------------------------
         % 2. Initialisation of variables to obtain in the postprocessing.
         %--------------------------------------------------------------------------
         switch str.postproc.output_type
             case 'averaging'
                 dim                                           =  str.data.dim;
                 nnode                                         =  str.n_nodes;
                 %-----------------------------------------
                 % Second Piola and the different contributions.
                 %-----------------------------------------
                 str.postproc.stress                           =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_mec_iso                   =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_mec_aniso                 =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_mec                       =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_diel_iso                  =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_diel_aniso                =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_diel                      =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_piezo                     =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_elec                      =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_deviatoric                =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.stress_volumetric                =  zeros((dim+1)*dim/2,nnode);
                 %-----------------------------------------
                 % Cauchy stress and the different contributions.
                 %-----------------------------------------
                 str.postproc.sigma                            =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_mec_iso                    =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_mec_aniso                  =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_mec                        =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_diel_iso                   =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_diel_aniso                 =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_diel                       =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_piezo                      =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_elec                       =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_deviatoric                 =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_volumetric                 =  zeros((dim+1)*dim/2,nnode);
                 str.postproc.sigma_pressure                   =  zeros(nnode,1);
                 %-----------------------------------------
                 % Eulerian and Lagrangian electric field.
                 %-----------------------------------------
                 str.postproc.E0                               =  zeros(dim,nnode);
                 str.postproc.E                                =  zeros(dim,nnode);
                 %-----------------------------------------
                 % Eulerian and Lagrangian electric displacement.
                 %-----------------------------------------
                 str.postproc.D                                =  zeros(dim,nnode);
                 str.postproc.D_diel_iso                       =  zeros(dim,nnode);
                 str.postproc.D_diel_aniso                     =  zeros(dim,nnode);
                 str.postproc.D_diel                           =  zeros(dim,nnode);
                 str.postproc.D_piezo                          =  zeros(dim,nnode);
                 str.postproc.D0                               =  zeros(dim,nnode);
                 str.postproc.D0_diel_iso                      =  zeros(dim,nnode);
                 str.postproc.D0_diel_aniso                    =  zeros(dim,nnode);
                 str.postproc.D0_diel                          =  zeros(dim,nnode);
                 str.postproc.D0_piezo                         =  zeros(dim,nnode);
                 %-----------------------------------------
                 % Deformation gradient and jacobian.
                 %-----------------------------------------
                 str.postproc.F                                =  zeros(dim,dim,nnode);
                 str.postproc.J                                =  zeros(nnode,1);
                 %-----------------------------------------
                 % Plastic information.
                 %-----------------------------------------
                 str.postproc.hardening_variable               =  zeros(nnode,1);
                 str.postproc.lambda_p                         =  zeros(dim,nnode);
                 %-----------------------------------------
                 % Energy density function.
                 %-----------------------------------------
                 str.postproc.Phi                              =  zeros(nnode,1);
                 str.postproc.Phi_mec_iso                      =  zeros(nnode,1);
                 str.postproc.Phi_mec_aniso                    =  zeros(nnode,1);
                 str.postproc.Phi_mec                          =  zeros(nnode,1);
                 str.postproc.Phi_diel_iso                     =  zeros(nnode,1);
                 str.postproc.Phi_diel_aniso                   =  zeros(nnode,1);
                 str.postproc.Phi_diel                         =  zeros(nnode,1);
                 str.postproc.Phi_piezo                        =  zeros(nnode,1);
                 str.postproc.Phi_elec                         =  zeros(nnode,1);
                 %-----------------------------------------
                 % Electric density vector.
                 %-----------------------------------------
                 str.postproc.J_electric                       =  zeros(dim,nnode);
                 str.postproc.J0_electric                      =  zeros(dim,nnode);
             case 'discontinuous'
         end
         switch str.data.analysis
             case 'dynamic'
                 if str.time_iteration>=2     
                     str.postproc.Dn1                          =  str.postproc.D0;       % Electric displacement in the previous time step.
                 else
                     str.postproc.Dn1                          =  zeros(dim,str.n_nodes);         % Electric displacement in the previous time step.
                 end
         end

         %--------------------------------------------------------------------------
         % 3. This function stores the different component of the stresses in vector
         %    form.
         %--------------------------------------------------------------------------
         f                                                     =  @(A,i,j) reshape(squeeze(A(i,j,:)),1,size(str.quadrature.Chi,1));
         switch str.data.dim
             case 2
                 g                                             =  @(B) [f(B,1,1);f(B,2,2);f(B,1,2)];
             case 3     
                 g                                             =  @(B) [f(B,1,1);f(B,2,2);f(B,3,3);f(B,1,2);f(B,1,3);f(B,2,3)];
         end

         %-----------------------------------------------------------
         % 4. Get the variables to plot in the nodes.
         %-----------------------------------------------------------
         str.nodes_counter                                     = zeros(1,str.n_nodes);
         for ielem=1:str.n_elem
             %nodes_elem_plot                                  =  str.connectivity_plot(ielem,:); This is the original version.
             nodes_elem_plot                                   =  str.connectivity(ielem,:);
             nodes_elem                                        =  str.connectivity(ielem,:);
             str.properties.material_identifier                =  str.properties.material(ielem,1);
             %-----------------------------
             % Current time step.
             %-----------------------------
             Xelem                                             =  str.Lagrangian_X(:,nodes_elem);
             xelem                                             =  str.Eulerian_x(:,nodes_elem);
             phielem                                           =  str.phi(nodes_elem);
             str                                               =  gradients(xelem,Xelem,phielem,str);
             Jn                                                =  str.grad.J;
             Fn                                                =  str.grad.F;
             %-----------------------------------------
             % Second Piola and different component.
             %-----------------------------------------
             if 0
             [S,S_mec_iso,S_mec_aniso,...
                 S_mec,S_diel_iso,S_diel_aniso,...
                 S_diel,S_piezo,S_elec]                        =  second_piola(str);
             end
             switch str.vacuum.flag
                    case 0 
                         switch str.formulation
                             case 'D_to_E'
                                  [S]                          =  second_piola_d(str);
                             case 'E'
                                  [S]                          =  second_piola(str);
                         end                                 
                    case 1
                         vacuum_identifier                     =  max(str.properties.material);
                         switch str.properties.material(ielem) 
                                case vacuum_identifier
                                     str.formulation           =  'E';
                                     [S]                       =  second_piola_vacuum(str);
                                otherwise
                                     str.formulation           =  str.old_formulation;
                                     [S]                       =  second_piola_d(str);
                         end
             end  
             %-----------------------------------------
             % Lagrangian electric displacement.
             %-----------------------------------------
             if 0
             [D0,D0_diel_iso,D0_diel_aniso,...     
                 D0_diel,D0_piezo]                             =  electric_displacement(str);
             end
             str.ielem                                         =  ielem;
             switch str.vacuum.flag
                    case 0
                         switch str.formulation
                             case 'D_to_E'
                                  [D0,str]                     =  internal_Newton_Raphson(str);
                             case 'E'
                                  D0                           =  electric_displacement(str);
                         end
                    case 1
                         vacuum_identifier                     =  max(str.properties.material);
                         switch str.properties.material(ielem) 
                                case vacuum_identifier
                                     str.formulation           =  'E';
                                     [D0]                      =  electric_displacement_vacuum(str);
                             otherwise
                                     str.formulation           =  str.old_formulation;                                    
                                     [D0,str]                  =  internal_Newton_Raphson(str);
                         end
             end             
             %-----------------------------------------
             % Cauchy stress and different component.
             %-----------------------------------------
             sigma                                             =  cauchy_stress(str,S);
             if 0
             sigma_mec_iso                                     =  cauchy_stress(str,S_mec_iso);
             sigma_mec_aniso                                   =  cauchy_stress(str,S_mec_aniso);
             sigma_mec                                         =  cauchy_stress(str,S_mec);
             sigma_diel_iso                                    =  cauchy_stress(str,S_diel_iso);
             sigma_diel_aniso                                  =  cauchy_stress(str,S_diel_aniso);
             sigma_diel                                        =  cauchy_stress(str,S_diel);
             sigma_piezo                                       =  cauchy_stress(str,S_piezo);
             sigma_elec                                        =  cauchy_stress(str,S_elec);
             end
             %-----------------------------------------
             % Volumetric and deviatoric decomposition
             % of Cauchy and Piola.
             %-----------------------------------------
             [sigma_deviatoric,sigma_volumetric,...
                 S_deviatoric,S_volumetric,p]                  =  volumetric_deviatoric_cauchy_Piola_decomposition(str,sigma);
             %-----------------------------------------
             % Eulerian electric displacement and different components.
             %-----------------------------------------
             [eul_D]                                           =  elec_disp(str,D0);
             if 0
             eul_D_diel_iso                                    =  elec_disp(str,D0_diel_iso);
             eul_D_diel_aniso                                  =  elec_disp(str,D0_diel_aniso);
             eul_D_diel                                        =  elec_disp(str,D0_diel);
             eul_D_piezo                                       =  elec_disp(str,D0_piezo);
             end
             %----------------------------------------------------------------
             % Second Piola and different contributions.
             %----------------------------------------------------------------
             str.postproc.stress(:,nodes_elem_plot)            =  str.postproc.stress(:,nodes_elem_plot) + g(S);
             if 0
             str.postproc.stress_mec_iso(:,nodes_elem_plot)    =  str.postproc.stress_mec_iso(:,nodes_elem_plot) + g(S_mec_iso);
             str.postproc.stress_mec_aniso(:,nodes_elem_plot)  =  str.postproc.stress_mec_aniso(:,nodes_elem_plot) + g(S_mec_aniso);
             str.postproc.stress_mec(:,nodes_elem_plot)        =  str.postproc.stress_mec(:,nodes_elem_plot) + g(S_mec);
             str.postproc.stress_diel_iso(:,nodes_elem_plot)   =  str.postproc.stress_diel_iso(:,nodes_elem_plot) + g(S_diel_iso);
             str.postproc.stress_diel_aniso(:,nodes_elem_plot) =  str.postproc.stress_diel_aniso(:,nodes_elem_plot) + g(S_diel_aniso);
             str.postproc.stress_diel(:,nodes_elem_plot)       =  str.postproc.stress_diel(:,nodes_elem_plot) + g(S_diel);
             str.postproc.stress_piezo(:,nodes_elem_plot)      =  str.postproc.stress_piezo(:,nodes_elem_plot) + g(S_piezo);
             str.postproc.stress_elec(:,nodes_elem_plot)       =  str.postproc.stress_elec(:,nodes_elem_plot) + g(S_elec);
             end
             str.postproc.stress_deviatoric(:,nodes_elem_plot) =  str.postproc.stress_deviatoric(:,nodes_elem_plot) + g(S_deviatoric);
             str.postproc.stress_volumetric(:,nodes_elem_plot) =  str.postproc.stress_volumetric(:,nodes_elem_plot) + g(S_volumetric);
             %----------------------------------------------------------------
             % Cauchy stress and different contributions.
             %----------------------------------------------------------------
             str.postproc.sigma(:,nodes_elem_plot)             =  str.postproc.sigma(:,nodes_elem_plot) +  g(sigma);
             if 0
             str.postproc.sigma_mec_iso(:,nodes_elem_plot)     =  str.postproc.sigma_mec_iso(:,nodes_elem_plot) + g(sigma_mec_iso);
             str.postproc.sigma_mec_aniso(:,nodes_elem_plot)   =  str.postproc.sigma_mec_aniso(:,nodes_elem_plot) + g(sigma_mec_aniso);
             str.postproc.sigma_mec(:,nodes_elem_plot)         =  str.postproc.sigma_mec(:,nodes_elem_plot) + g(sigma_mec);
             str.postproc.sigma_diel_iso(:,nodes_elem_plot)    =  str.postproc.sigma_diel_iso(:,nodes_elem_plot) + g(sigma_diel_iso);
             str.postproc.sigma_diel_aniso(:,nodes_elem_plot)  =  str.postproc.sigma_diel_aniso(:,nodes_elem_plot) + g(sigma_diel_aniso);
             str.postproc.sigma_diel(:,nodes_elem_plot)        =  str.postproc.sigma_diel(:,nodes_elem_plot) + g(sigma_diel);
             str.postproc.sigma_piezo(:,nodes_elem_plot)       =  str.postproc.sigma_piezo(:,nodes_elem_plot) + g(sigma_piezo);
             str.postproc.sigma_elec(:,nodes_elem_plot)        =  str.postproc.sigma_elec(:,nodes_elem_plot) + g(sigma_elec);
             str.postproc.sigma_deviatoric(:,nodes_elem_plot)  =  str.postproc.sigma_deviatoric(:,nodes_elem_plot) + g(sigma_deviatoric);
             end
             str.postproc.sigma_volumetric(:,nodes_elem_plot)  =  str.postproc.sigma_volumetric(:,nodes_elem_plot) + g(sigma_volumetric);
             str.postproc.sigma_pressure(nodes_elem_plot,1)    =  str.postproc.sigma_pressure(nodes_elem_plot,1) + p;
             %----------------------------------------------------------------
             % Electric displacement. Eulerian part.
             %----------------------------------------------------------------
             str.postproc.D(:,nodes_elem_plot)                 =  str.postproc.D(:,nodes_elem_plot) + eul_D;
             if 0
             str.postproc.D_diel_iso(:,nodes_elem_plot)        =  str.postproc.D_diel_iso(:,nodes_elem_plot) + eul_D_diel_iso;
             str.postproc.D_diel_aniso(:,nodes_elem_plot)      =  str.postproc.D_diel_aniso(:,nodes_elem_plot) + eul_D_diel_aniso;
             str.postproc.D_diel(:,nodes_elem_plot)            =  str.postproc.D_diel(:,nodes_elem_plot) + eul_D_diel;
             str.postproc.D_piezo(:,nodes_elem_plot)           =  str.postproc.D_piezo(:,nodes_elem_plot) + eul_D_piezo;
             end
             %----------------------------------------------------------------
             % Electric displacement. Lagrangian part.
             %----------------------------------------------------------------
             str.postproc.D0(:,nodes_elem_plot)                =  str.postproc.D0(:,nodes_elem_plot) + D0;
             if 0
             str.postproc.D0_diel_iso(:,nodes_elem_plot)       =  str.postproc.D0_diel_iso(:,nodes_elem_plot) + D0_diel_iso;
             str.postproc.D0_diel_aniso(:,nodes_elem_plot)     =  str.postproc.D0_diel_aniso(:,nodes_elem_plot) + D0_diel_aniso;
             str.postproc.D0_diel(:,nodes_elem_plot)           =  str.postproc.D0_diel(:,nodes_elem_plot) + D0_diel;
             str.postproc.D0_piezo(:,nodes_elem_plot)          =  str.postproc.D0_piezo(:,nodes_elem_plot) + D0_piezo;
             end
             %----------------------------------------------------------------
             % Electric field.
             %----------------------------------------------------------------
             str.postproc.E0(:,nodes_elem_plot)                =  str.postproc.E0(:,nodes_elem_plot) + str.grad.E0;
             str.postproc.E(:,nodes_elem_plot)                 =  str.postproc.E(:,nodes_elem_plot) + str.grad.E;
             %----------------------------------------------------------------
             % Deformation.
             %----------------------------------------------------------------
             str.postproc.F(:,:,nodes_elem_plot)               =  str.postproc.F(:,:,nodes_elem_plot) + str.grad.F;
             str.postproc.J(nodes_elem_plot,1)                 =  str.postproc.J(nodes_elem_plot,1) + str.grad.J;
             %----------------------------------------------------------------
             % Stored energy density function.
             %----------------------------------------------------------------
             if 0
             [str]                                             =  stored_energy_density_function(str);
             str.postproc.Phi(nodes_elem_plot,1)               =  str.postproc.Phi(nodes_elem_plot,1) + str.Phi;
             str.postproc.Phi_mec_iso(nodes_elem_plot,1)       =  str.postproc.Phi_mec_iso(nodes_elem_plot,1) + str.Phi_mec_iso;
             str.postproc.Phi_mec_aniso(nodes_elem_plot,1)     =  str.postproc.Phi_mec_aniso(nodes_elem_plot,1) + str.Phi_mec_aniso;
             str.postproc.Phi_mec(nodes_elem_plot,1)           =  str.postproc.Phi_mec(nodes_elem_plot,1) + str.Phi_mec;
             str.postproc.Phi_diel_iso(nodes_elem_plot,1)      =  str.postproc.Phi_diel_iso(nodes_elem_plot,1) + str.Phi_diel_iso;
             str.postproc.Phi_diel_aniso(nodes_elem_plot,1)    =  str.postproc.Phi_diel_aniso(nodes_elem_plot,1) + str.Phi_diel_aniso;
             str.postproc.Phi_diel(nodes_elem_plot,1)          =  str.postproc.Phi_diel(nodes_elem_plot,1) + str.Phi_diel;
             str.postproc.Phi_piezo(nodes_elem_plot,1)         =  str.postproc.Phi_piezo(nodes_elem_plot,1) + str.Phi_piezo;
             str.postproc.Phi_elec(nodes_elem_plot,1)          =  str.postproc.Phi_elec(nodes_elem_plot,1) + str.Phi_elec;
             end
             %-----------------------------
             % Previous time step in order to get J.
             %-----------------------------
             switch str.data.analysis
                 case 'dynamic'
                     if str.time_iteration>=2
                         switch str.formulation
                             case 'E'
                                  xelemn1                      =  str.Eulerian_xn1(:,nodes_elem);
                                  phielemn1                    =  str.phin1(nodes_elem);
                                  [str]                        =  gradients(xelemn1,Xelem,phielemn1,str);
                                  [D0n1]                       =  electric_displacement(str);
                                  [J_electric,J0_electric]     =  Intensity_vector(str,D0,D0n1,Jn,Fn);
                                  str.postproc.J_electric(:,...
                                      nodes_elem_plot)         =  str.postproc.J_electric(:,nodes_elem_plot) + J_electric;
                                  str.postproc.J0_electric(:,...
                                      nodes_elem_plot)         =  str.postproc.J0_electric(:,nodes_elem_plot) + J0_electric;
                             case 'D_to_E'
                                  xelemn1                      =  str.Eulerian_xn1(:,nodes_elem);
                                  phielemn1                    =  str.phin1(nodes_elem);
                                  [str]                        =  gradients(xelemn1,Xelem,phielemn1,str);
%                                 [D0n1]                       =  electric_displacement(str);                                  
                                  [D0n1,str]                   =  internal_Newton_Raphson(str);                                  
                                  [J_electric,J0_electric]     =  Intensity_vector(str,D0,D0n1,Jn,Fn);
                                  str.postproc.J_electric(:,...
                                      nodes_elem_plot)         =  str.postproc.J_electric(:,nodes_elem_plot) + J_electric;
                                  str.postproc.J0_electric(:,...
                                      nodes_elem_plot)         =  str.postproc.J0_electric(:,nodes_elem_plot) + J0_electric;                                 
                         end
                     end
             end
             switch str.data.elastoplasticity
                 case 1
                      [S,Dv,new_str]                                       =  Radial_Return_Mapping_algorithm(str);
                      for inode=1:size(str.connectivity,2)
                          lambda_p                                         =  eig(str.solid.kinematic.plastic.C_p(:,:,ielem,inode));
                          lambda_p                                         =  sqrt(lambda_p);
                          node                                             =  str.connectivity(ielem,inode);
                          str.postproc.hardening_variable(node)            =  str.postproc.hardening_variable(inode) + new_str.solid.kinematic.plastic.hardening_variable(ielem,inode);
                          str.postproc.lambda_p(:,node)                    =  str.postproc.lambda_p(:,node)  +  lambda_p;
                      end
                 case 0
                      str.ep(nodes_elem_plot)                  =  str.ep(nodes_elem_plot) + zeros(dim,dim,nnode);
                      str.Cp(nodes_elem_plot)                  =  str.Cp(nodes_elem_plot) + zeros(dim,dim,nnode);
                      for inode=1:size(str.connectivity,2)
                          lambda_p                                         =  sqrt(eig(eye(dim)));
                          str.postproc.lambda_p(:,nodes_elem_plot)         =  str.postproc.lambda_p(:,nodes_elem_plot)  +  lambda_p;
                      end
             end
             
             %----------------------------------------------------------------
             % Element contributions on each node.
             %----------------------------------------------------------------
             str.nodes_counter(nodes_elem_plot)                =  str.nodes_counter(nodes_elem_plot) + 1;
         end
         %--------------------------------------------------------------------------
         % 5. Remove the inner nodes from the variables.
         %--------------------------------------------------------------------------
         [str]                                                 =  plot_variables_treatment(str);

         %--------------------------------------------------------------------------
         % 6. Method to compute the variables in the vertices(output type): averaging...
         %--------------------------------------------------------------------------
         [str]                                                 =  output_type(str);
 
    case 1
         %--------------------------------------------------------------------------
         % 5. Remove the inner nodes from the variables.
         %--------------------------------------------------------------------------
         [str]                                                 =  plot_variables_treatment(str);
         %--------------------------------------------------
         % This case is for when we need to compare examples. 
         %--------------------------------------------------
         str                                                   =  comparison_postprocessing(str);        
end

% str.postproc.nodes = str.nodes';
% str.postproc.t = str.connectivity';

str.data.dim                                                    =  old_str.data.dim;