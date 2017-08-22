%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the postpocessing per element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  str                                                              =  postpocessing_computing_paraview(str)

str.f_e                                                                    =  str.postproc.f_e;
str.quadrature.Chi                                                         =  str.nodal_positions_postprocessing_mesh;
%--------------------------------------------------------------------------
% This function stores the different component of the stresses in vector
% form.
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
%--------------------------------------------------------------------------
% Get the variables to plot in the nodes. 
%--------------------------------------------------------------------------
for ielem=1:str.n_elem
    %fprintf('%d\n',ielem)
    str.ielem                                                              =  ielem;
    nodes_elem_plot                                                        =  str.postproc.connectivity(ielem,:);
    solution_nodes_elem                                                    =  str.connectivity(ielem,:);
    str.properties.material_identifier                                     =  str.properties.material(ielem,1);
    %----------------------------------------------------------------------
    % Current time step.
    %----------------------------------------------------------------------
    gauss_level_information                                                =  postprocessing_element_nodes_paraview(str);    
    Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem);        
    xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem);
    phielem                                                                =  str.phi(solution_nodes_elem,1);        
    
    xelem_postproc                                                         =  gauss_level_information.Eulerian_x;
    Xelem_postproc                                                         =  gauss_level_information.Lagrangian_X;
    phielem_postproc                                                       =  gauss_level_information.phi;    
    str                                                                    =  gradients(xelem,Xelem,phielem,str);
    %----------------------------------------------------------------------
    % x and phi postprocessing nodes. 
    %----------------------------------------------------------------------
    str.postproc.Eulerian_x(:,nodes_elem_plot)                             =  xelem_postproc;
    str.postproc.displacement(:,nodes_elem_plot)                           =  xelem_postproc - Xelem_postproc;
    str.postproc.phi(nodes_elem_plot)                                      =  phielem_postproc';
    %----------------------------------------------------------------------
    % Initialising everything to zero
    %----------------------------------------------------------------------
    F                                                                      =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    H                                                                      =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    J                                                                      =  zeros(size(str.postproc.f_e.N,2),1);
    SigmaF                                                                 =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    SigmaH                                                                 =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    SigmaJ                                                                 =  zeros(size(str.postproc.f_e.N,2),1);
    E0                                                                     =  zeros(dim,size(str.postproc.f_e.N,2));
    D0                                                                     =  zeros(dim,size(str.postproc.f_e.N,2));
    Kirchhoff                                                              =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
    p                                                                      =  zeros(size(str.postproc.f_e.N,2),1);
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
            [F,H,J,E0,D0,Kirchhoff]        =  displacement_potential_postprocessing(str,ielem);
        %------------------------------------------------------------------
        %  Mixed formulations  
        %------------------------------------------------------------------    
        case 'mixed_formulation'            
             str.properties.stabilisation_parameter_momentum=0;
             gauss_level_information                                       =  gauss_level_information_mixed_formulations(str);
             %gauss_level_information                                      =  gauss_level_information_mixed_formulations_final(str);
             switch str.data.formulation.mixed_type 
                    case 'displacement_formulation' 
                         [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff]            =  displacement_postprocessing(SigmaF,SigmaH,SigmaJ,Kirchhoff,gauss_level_information);
                    case 'full_mixed_formulation' 
                         [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff]            =  F_H_J_SF_SH_SJ_postprocessing(gauss_level_information,str);                        
                    case {'complementary_energy','stabilised_complementary_energy','stabilised_complementary_energy_1','stabilised_complementary_energy_2','stabilised_complementary_energy_3'}
                         % Not defined
                    case {'full_incompressible','compressible_two_field'}
                         [F,H,J,Kirchhoff,p_formulation]                   =  two_field_mechanical_postprocessing(str);
                    case {'superfull_mixed_incompressible_formulation','full_mixed_incompressible_formulation'}
                        % Not defined
                    case 'displacement_potential_mixed_formulation'
                         [F,H,J,E0,D0,SigmaF,SigmaH,SigmaJ,Kirchhoff]      =  displacement_potential_current_postprocessing(str,gauss_level_information,Kirchhoff,...
                                                                                 SigmaF,SigmaH,SigmaJ,D0);                        
                    case 'u_p_phi_incompressible'  
                         [F,H,J,E0,D0,SigmaF,SigmaH,SigmaJ,Kirchhoff,p]    =  u_p_phi_incompressible_postprocessing(str,gauss_level_information,SigmaF,SigmaH,SigmaJ,Kirchhoff,D0,p);
                         p_formulation                                     =  p;
                    case 'full_mixed_formulation_electroelasticity_F_D0'
                         [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff,E0,D0]      =  F_D0_postprocessing(str,gauss_level_information);
                    case 'full_mixed_formulation_electroelasticity_F_D0_V'  
                         [option,F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff,...
                             D0,E0]                                        =  F_D0_V_postprocessing(str,gauss_level_information,E0);    
                    case 'full_mixed_formulation_electroelasticity_F_D0_V_incompressible'  
                         [option,F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff,...
                             D0,E0,p_formulation]                          =  F_D0_V_p_postprocessing(str,gauss_level_information,E0);    
             end
    end
    %----------------------------------------------------------------------            
    % Cauchy stress and different component.
    %----------------------------------------------------------------------            
    [sigma]                                                              =  cauchy_stress_paraview(str,Kirchhoff,F,H,J,SigmaF,SigmaH,SigmaJ,option,gauss_level_information);
    %[sigma,p]                                                              =  cauchy_stress_paraview(str,Kirchhoff,F,H,J,SigmaF,SigmaH,SigmaJ,option,gauss_level_information);
    %----------------------------------------------------------------------            
    % Update pressure
    %----------------------------------------------------------------------            
    switch str.data.formulation.mixed_type
        case {'full_incompressible','compressible_two_field','nearly_incompressible_three_field','nearly_incompressible_two_field',...
                'full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation','u_p_phi_compressible',...
                'u_p_phi_incompressible','full_mixed_formulation_electroelasticity_F_D0_V_incompressible'}            
             %p                                                             =  p + p_formulation;
    end    
    %----------------------------------------------------------------------            
    % Eulerian electric displacement and different components.
    %----------------------------------------------------------------------
    [eul_D,E]                                                              =  elec_disp_paraview(str,D0,E0);
    %----------------------------------------------------------------------            
    % First Piola stress conjugates. 
    %----------------------------------------------------------------------            
    str.postproc.SigmaF(:,:,nodes_elem_plot)                               =  str.postproc.SigmaF(:,:,nodes_elem_plot) + twoD_treatment(SigmaF);
    %str.SigmaFF(:,ielem)  =  reshape(SigmaF(1,1,:),[],1);
    str.postproc.SigmaH(:,:,nodes_elem_plot)                               =  str.postproc.SigmaH(:,:,nodes_elem_plot) + twoD_treatment(SigmaH);
    str.postproc.SigmaJ(nodes_elem_plot,1)                                 =  str.postproc.SigmaJ(nodes_elem_plot,1) + SigmaJ;
    str.postproc.First_Piola(:,:,nodes_elem_plot)                          =  str.postproc.First_Piola(:,:,nodes_elem_plot) + twoD_treatment(Kirchhoff);
    %----------------------------------------------------------------------            
    % Cauchy stress and different contributions.
    %----------------------------------------------------------------------            
    str.postproc.sigma(:,nodes_elem_plot)                                  =  str.postproc.sigma(:,nodes_elem_plot) +  g(sigma);
    str.postproc.sigma_pressure(nodes_elem_plot,1)                         =  str.postproc.sigma_pressure(nodes_elem_plot,1) + p;
    %----------------------------------------------------------------------            
    % Electric displacement. Eulerian part.
    %----------------------------------------------------------------------            
    str.postproc.D(:,nodes_elem_plot)                                      =  str.postproc.D(:,nodes_elem_plot) + eul_D;
    %----------------------------------------------------------------------            
    % Electric displacement. Lagrangian part.
    %----------------------------------------------------------------------            
    str.postproc.D0(:,nodes_elem_plot)                                     =  str.postproc.D0(:,nodes_elem_plot) + D0;
    %----------------------------------------------------------------------            
    % Electric field.
    %----------------------------------------------------------------------            
    str.postproc.E0(:,nodes_elem_plot)                                     =  str.postproc.E0(:,nodes_elem_plot) + E0;
    str.postproc.E(:,nodes_elem_plot)                                      =  str.postproc.E(:,nodes_elem_plot) + E;
    %----------------------------------------------------------------------            
    % Deformation.
    %----------------------------------------------------------------------            
    str.postproc.F(:,:,nodes_elem_plot)                                    =  str.postproc.F(:,:,nodes_elem_plot) + twoD_treatment(F);
    str.postproc.H(:,:,nodes_elem_plot)                                    =  str.postproc.H(:,:,nodes_elem_plot) + twoD_treatment(H);
    str.postproc.J(nodes_elem_plot,1)                                      =  str.postproc.J(nodes_elem_plot,1) + J;
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

