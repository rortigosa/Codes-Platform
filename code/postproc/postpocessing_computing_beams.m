%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the postpocessing per element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  str                                                              =  postpocessing_computing_beams(str)

str.f_e                                                                    =  str.solid.BEAM_SHELL.continuum.isoparametric_f_e;
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
    %nodes_elem_plot                                                       =  str.postproc.connectivity(ielem,:);
    solution_nodes_elem                                                    =  str.connectivity(ielem,:);
    str.properties.material_identifier                                     =  str.properties.material(ielem,1);
    %-----------------------------
    % Current time step.
    %-----------------------------
    gauss_level_information                                                =  postprocessing_element_nodes(str);   
    Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem) + 0*repmat([2;2;0],1,8);        
    xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem) + 0*repmat([2;2;0],1,8);
%     if ielem<=160
%     Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem);        
%     xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem);
%     end
%     if ielem>160
%     Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem) + 0*repmat([0;0;-0.1],1,8);        
%     xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem) + 0*repmat([0;0;-0.1],1,8);
%     end
%     if ielem>320
%     Xelem                                                                  =  str.Lagrangian_X(:,solution_nodes_elem) + 0*repmat([0;0;-0.2],1,8);        
%     xelem                                                                  =  str.Eulerian_x(:,solution_nodes_elem) + 0*repmat([0;0;-0.2],1,8);
%     end    
    phielem                                                                =  str.phi(solution_nodes_elem,1);        
    velocity                                                               =  reshape(str.velocity,str.data.dim,str.n_nodes);
    velocity                                                               =  velocity(:,str.connectivity(ielem,:));
    acceleration                                                           =  reshape(str.acceleration,str.data.dim,str.n_nodes);
    acceleration                                                           =  acceleration(:,str.connectivity(ielem,:));
    str                                                                    =  gradients(xelem,Xelem,phielem,str);
    %-----------------------------
    % x and phi postprocessing nodes.
    %-----------------------------        
    str.postproc.Eulerian_x(:,solution_nodes_elem)                         =  xelem;
    str.postproc.phi(solution_nodes_elem)                                  =  phielem';
    str.postproc.velocity(:,solution_nodes_elem)                           =  velocity;
    str.postproc.acceleration(:,solution_nodes_elem)                       =  acceleration;
    %----------------------------------------
    % Initialising everything to zero
    %----------------------------------------
    F                                                                      =  zeros(dim,dim,size(str.f_e.N,2));
    H                                                                      =  zeros(dim,dim,size(str.f_e.N,2));
    J                                                                      =  zeros(size(str.f_e.N,2),1);
    SigmaF                                                                 =  zeros(3,3,size(str.f_e.N,2));
    SigmaH                                                                 =  zeros(3,3,size(str.f_e.N,2));
    SigmaJ                                                                 =  zeros(size(str.f_e.N,2),1);
    E0                                                                     =  zeros(dim,size(str.f_e.N,2));
    E                                                                      =  zeros(dim,size(str.f_e.N,2));
    D0                                                                     =  zeros(dim,size(str.f_e.N,2));
    S                                                                      =  zeros(dim,dim,size(str.f_e.N,2));    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Computation of variables for the vacuum and solid meshes (Only in BEM)
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    switch str.data.formulation_type          
        %------------------------------------------------------------------
        %  Displacement potential formulation
        %------------------------------------------------------------------
        case 'displacement_potential_formulation' 
             oldstr                                                        =  str;
             new_str                                                       =  gradients(xelem,Xelem,phielem,str);
             str                                                           =  oldstr;
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
             switch str.formulation
                 case 'D_to_E'
                      [D0,str]                                             =  internal_Newton_Raphson(str);
                 case 'E'
                      D0                                                   =  electric_displacement(str);
             end
    end
    %----------------------------------------------------------------------
    % Compute norms of tensors
    %----------------------------------------------------------------------    
    FL2norm                                                                =  zeros(size(str.f_e.N,2),1);
    for igauss=1:size(str.f_e.N,2)
        FL2norm(igauss)                                                    =  trace(F(:,:,igauss)'*F(:,:,igauss));
    end
    
    %----------------------------------------------------------------------
    % Eigenvalues of the second derivative of the energy functional with
    % respect to the deformation gradient F. This is done for stability
    % analysis purposes.
    %---------------------------------------------------------------------- 
    eigenvalues                                                            =  zeros(str.data.dim^2,size(str.quadrature.Chi,1));
    %-----------------------------------------
    % Cauchy stress and different component.
    %-----------------------------------------
    sigma                                                                  =  cauchy_stress(str,Kirchhoff,F,H,J,SigmaF,SigmaH,SigmaJ,gauss_level_information);
    %-----------------------------------------
    % Cauchy stress and different component.
    %-----------------------------------------
    [SigmaF,SigmaH,SigmaJ]                                                 =  conjugate_stresses(str,F,H,J);
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
                'full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation'}
            
             p                                                             =  p_old';
    end    
    %-----------------------------------------            
    % Eulerian electric displacement and different components.
    %-----------------------------------------
    [eul_D]                                                                =  elec_disp(str,D0);
    %----------------------------------------------------------------
    % Second Piola and different contributions. 
    %----------------------------------------------------------------
    str.postproc.stress(:,solution_nodes_elem)                             =  str.postproc.stress(:,solution_nodes_elem) + g(S);
    str.postproc.stress_deviatoric(:,solution_nodes_elem)                  =  str.postproc.stress_deviatoric(:,solution_nodes_elem) + g(S_deviatoric);
    str.postproc.stress_volumetric(:,solution_nodes_elem)                  =  str.postproc.stress_volumetric(:,solution_nodes_elem) + g(S_volumetric);
    %----------------------------------------------------------------
    % First Piola stress conjugates.
    %----------------------------------------------------------------
    str.postproc.SigmaF(:,:,solution_nodes_elem)                           =  str.postproc.SigmaF(:,:,solution_nodes_elem) + twoD_treatment(SigmaF);
    str.postproc.SigmaH(:,:,solution_nodes_elem)                           =  str.postproc.SigmaH(:,:,solution_nodes_elem) + twoD_treatment(SigmaH);
    str.postproc.SigmaJ(solution_nodes_elem,1)                             =  str.postproc.SigmaJ(solution_nodes_elem,1) + SigmaJ;
    str.postproc.First_Piola(:,:,solution_nodes_elem)                      =  str.postproc.First_Piola(:,:,solution_nodes_elem) + twoD_treatment(Kirchhoff);
    %----------------------------------------------------------------
    % Cauchy stress and different contributions.
    %----------------------------------------------------------------
    str.postproc.sigma(:,solution_nodes_elem)                              =  str.postproc.sigma(:,solution_nodes_elem) +  g(sigma);
    str.postproc.sigma_deviatoric(:,solution_nodes_elem)                   =  str.postproc.sigma_deviatoric(:,solution_nodes_elem) + g(sigma_deviatoric);
    str.postproc.sigma_volumetric(:,solution_nodes_elem)                   =  str.postproc.sigma_volumetric(:,solution_nodes_elem) + g(sigma_volumetric);
    str.postproc.sigma_pressure(solution_nodes_elem,1)                     =  str.postproc.sigma_pressure(solution_nodes_elem,1) + p;
    %----------------------------------------------------------------
    % Electric displacement. Eulerian part.
    %----------------------------------------------------------------
    str.postproc.D(:,solution_nodes_elem)                                  =  str.postproc.D(:,solution_nodes_elem) + eul_D;
    %----------------------------------------------------------------
    % Electric displacement. Lagrangian part.
    %----------------------------------------------------------------
    str.postproc.D0(:,solution_nodes_elem)                                 =  str.postproc.D0(:,solution_nodes_elem) + D0;
    %----------------------------------------------------------------
    % Electric field.
    %----------------------------------------------------------------
    str.postproc.E0(:,solution_nodes_elem)                                 =  str.postproc.E0(:,solution_nodes_elem) + E0;
    str.postproc.E(:,solution_nodes_elem)                                  =  str.postproc.E(:,solution_nodes_elem) + E;
    %----------------------------------------------------------------
    % Deformation.
    %----------------------------------------------------------------
    str.postproc.F(:,:,solution_nodes_elem)                                =  str.postproc.F(:,:,solution_nodes_elem) + twoD_treatment(F);
    str.postproc.FL2norm(solution_nodes_elem)                              =  str.postproc.FL2norm(solution_nodes_elem) + FL2norm;                  
    str.postproc.H(:,:,solution_nodes_elem)                                =  str.postproc.H(:,:,solution_nodes_elem) + twoD_treatment(H);
    str.postproc.J(solution_nodes_elem,1)                                  =  str.postproc.J(solution_nodes_elem,1) + J;
    %----------------------------------------------------------------
    % Second derivative of the energy function with respect to F. 
    %----------------------------------------------------------------
    str.postproc.eigenvalues(:,solution_nodes_elem)                        =  str.postproc.eigenvalues(:,solution_nodes_elem) + eigenvalues;
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