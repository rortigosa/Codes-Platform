function   [Kinetic_energy,Elastic_energy,Total_energy,...
    Kinetic_power,Elastic_power,Total_power]        =  conservation_energy_monitoring(str)

Kinetic_energy                                      =  0;
Elastic_energy                                      =  0;
Kinetic_power                                       =  0;
Elastic_power                                       =  0;

velocity                                            =  reshape(str.velocity,str.data.dim,size(str.velocity,1)/str.data.dim);
acceleration                                        =  reshape(str.acceleration,str.data.dim,size(str.acceleration,1)/str.data.dim);
for ielem=1:str.n_elem
    str.ielem                                       =  ielem;
    %----------------------------------------------------------------------
    % Compute information at gauss level.
    %----------------------------------------------------------------------
    nodes                                           =  str.connectivity(ielem,:);
    xelem                                           =  str.Eulerian_x(:,nodes);
    Xelem                                           =  str.Lagrangian_X(:,nodes);
    phielem                                         =  str.phi(nodes,1);
    str                                             =  gradients(xelem,Xelem,phielem,str);
    velocity_element                                =  velocity(:,nodes);
    velocity_gauss                                  =  velocity_element*str.f_e.N;
    acceleration_element                            =  acceleration(:,nodes);
    acceleration_gauss                              =  acceleration_element*str.f_e.N;
    gauss_level_information                         =  gauss_level_information_mixed_formulations(str);
    Piola                                           =  First_Piola_Kirchhoff_stress_tensor(gauss_level_information,str);    
    switch str.data.formulation.mixed_type
        case {'displacement_formulation','full_incompressible','full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation','compressible_two_field'}
            F                                       =  str.grad.F;
            H                                       =  str.grad.H;
            J                                       =  str.grad.J;
        case {'full_mixed_formulation'}
            F                                       =  gauss_level_information.F;
            H                                       =  gauss_level_information.H;
            J                                       =  gauss_level_information.J;
        case 'H_J_mixed_formulation'
            F                                       =  str.grad.F;
            H                                       =  gauss_level_information.H;
            J                                       =  gauss_level_information.J;
        case 'complementary_energy'        
            F                                       =  str.grad.F;  %  It't not true but it's ok momentaneously
            H                                       =  str.grad.H;
            J                                       =  str.grad.J;
    end
    for igauss=1:size(str.quadrature.Chi,1) 
        v_gauss                                     =  velocity_gauss(:,igauss);
        a_gauss                                     =  acceleration_gauss(:,igauss);
        F_3D                                        =  F(:,:,igauss);
        H_3D                                        =  H(:,:,igauss);
        J_3D                                        =  J(igauss);
        %------------------------------------------------------------------
        % Isoparametric information
        %------------------------------------------------------------------
        J_t                                         =  abs(det(str.grad.DX_chi(:,:,igauss)));
        W                                           =  str.quadrature.W_v(igauss);
        %------------------------------------------------------------------
        % Elastic energy
        %------------------------------------------------------------------
        switch str.material_model_info.material_model
            %--------------------------------------------------------------
            % Compressible Mooney-Rivlin in Javier's paper
            %--------------------------------------------------------------
            case 'simple_Mooney_Rivlin'
                 switch str.data.formulation.mixed_type
                        case {'full_incompressible'}
                             alpha                  =  str.properties.material_parameters.eta;
                             beta                   =  str.properties.material_parameters.gamma;
                             W_elastic              =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D)^(3/2) + gauss_level_information.p(igauss)*(J_3D-1);
                        case 'compressible_two_field'
                             alpha                  =  str.properties.material_parameters.eta;
                             beta                   =  str.properties.material_parameters.gamma;
                             lambda                 =  str.properties.material_parameters.k;
                             pressure               =  gauss_level_information.p(igauss);
                             %W_elastic              =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D)^(3/2) + lambda/2*(J_3D-1)^2 - 3*alpha - 3^(3/2)*beta;
                             W_elastic              =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D)^(3/2) + pressure*(J_3D-1) - 3*alpha - 3^(3/2)*beta;
                     otherwise
                             alpha                  =  str.properties.material_parameters.eta;
                             beta                   =  str.properties.material_parameters.gamma;
                             lambda                 =  str.properties.material_parameters.k;
                             %W_elastic             =  alpha*trace(F_3D'*F_3D) +  beta*trace(H_3D'*H_3D) - 2*(alpha+2*beta)*log(J_3D) + lambda*(J_3D-1)^2 - 3*(alpha+beta);                                                     
                             W_elastic              =  alpha*trace(F_3D'*F_3D) +  beta*trace(H_3D'*H_3D) - 4*beta*J_3D - 2*alpha*log(J_3D) + lambda*(J_3D-1)^2 - 3*(alpha+beta) + 4*beta;                                                     

                 end
            %--------------------------------------------------------------
            % Compressible Mooney-Rivlin in Javier's paper
            %--------------------------------------------------------------
            case 'simple_Mooney_Rivlin_v2'
                 switch str.data.formulation.mixed_type
                        case 'full_incompressible'
                             alpha                  =  str.properties.material_parameters.eta;
                             beta                   =  str.properties.material_parameters.gamma;
                             W_elastic              =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D) + gauss_level_information.p(igauss)*(J_3D-1);
                     otherwise
                             alpha                  =  str.properties.material_parameters.alpha;
                             beta                   =  str.properties.material_parameters.beta;
                             lambda                 =  str.properties.material_parameters.lambda;
                             W_elastic             =  alpha*trace(F_3D'*F_3D) +  beta*trace(H_3D'*H_3D) - 2*(alpha+2*beta)*log(J_3D) + lambda*(J_3D-1)^2 - 3*(alpha+beta);                                                     

                 end
                 %---------------------------------------------------------
                 % Compressible Mooney-Rivlin in Schroder's paper
                 %---------------------------------------------------------
            case 'Schroder_Mooney_Rivlin'
                 alpha                              =  str.properties.material_parameters.alpha;
                 beta                               =  str.properties.material_parameters.beta;
                 lambda                             =  str.properties.material_parameters.lambda;
                 epsilon                            =  str.properties.material_parameters.epsilon;
                 W_elastic                          =  alpha*trace(F_3D'*F_3D)^2 +  beta*trace(H_3D'*H_3D)^2 - 12*(alpha+2*beta)*log(J_3D) + ...
                                                       lambda/(2*epsilon^2)*(J_3D^epsilon + J_3D^(-epsilon) - 1) - 9*(alpha+beta);
                 %---------------------------------------------------------
                 % Deviatoric Mooney-Rivlin
                 %---------------------------------------------------------
            case 'Deviatoric_Mooney_Rivlin'
                 alpha                              =  str.properties.material_parameters.alpha;
                 beta                               =  str.properties.material_parameters.beta;
                 switch str.data.formulation.mixed_type
                     case {'full_mixed_incompressible_formulation','superfull_mixed_incompressible_formulation'}
                          W_elastic                          =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D)^(3/2) + gauss_level_information.p(igauss)*(J_3D-1) - 3*alpha - 3^(3/2)*beta;
                     otherwise
                          lambda                             =  str.properties.material_parameters.lambda;
                          W_elastic                          =  alpha*J_3D^(-2/3)*trace(F_3D'*F_3D) +  beta*J_3D^(-2)*trace(H_3D'*H_3D)^(3/2) + lambda*(J_3D-1)^2 - 3*alpha - 3^(3/2)*beta;
                 end                          
        end
        Elastic_energy                              =  Elastic_energy + W_elastic*W*J_t;
        DvDX                                        =  velocity_element*(str.grad.DN_X(:,:,igauss))';
        Elastic_power                               =  Elastic_power  + trace(Piola(:,:,igauss)'*DvDX)*(W*J_t);
        %------------------------------------------------------------------
        % Kinetic energy
        %------------------------------------------------------------------
        Kinetic_energy                              =  Kinetic_energy + 0.5*str.data.Rho*(v_gauss'*v_gauss)*W*J_t;
        Kinetic_power                               =  Kinetic_power  + str.data.Rho*(v_gauss'*a_gauss)*W*J_t;
    end
end
%------------------------------------------------------------------
% Total energy
%------------------------------------------------------------------
Total_energy                                        =  Elastic_energy + Kinetic_energy;
Total_power                                         =  Elastic_power  + Kinetic_power;
 