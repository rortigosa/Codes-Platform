
function minimum_q            =  acoustic_tensor_simulation(F,H,J,D0,d,n_alpha,n_beta,str)

plotting                      =  1;

q                             =  zeros(n_beta*n_alpha,1);
Time_N                        =  1;


switch plotting
    case 1
         [n_alpha,n_beta,Nodes_alpha,... 
          connectivity_alpha] =  spherical_parametrisation_mesh;
          data.control.outputfilename    =  fullfile('characterisation_experiments','results');
      
end  

 
for beta=1:n_beta
    for alpha=1:n_alpha
        theta_alpha           =  2*pi/n_alpha*(alpha-1);
        theta_beta            =  pi/n_beta*(beta-1);
        %------------------------------------------------------
        % Direction of propagation
        %------------------------------------------------------
        N                     =  [sin(theta_beta)*cos(theta_alpha);sin(theta_beta)*sin(theta_alpha);cos(theta_beta)];
        %------------------------------------------------------
        % Constitutive tensors
        %------------------------------------------------------
        [Conjugate,Hessian]   =  constitutive_law_conjugate_hessian_electro_acoustics(str.properties.material_parameters,F,H,J,D0,d,str.material_model_info.material_model);
        [C,QT,theta]          =  constitutive_tensors_electro_mechanics(F,H,J,D0,d,Conjugate,Hessian);
        %--------------------------------------------------------------
        % Acoustic tensor
        %--------------------------------------------------------------
        q(Time_N,1)           =  electro_acoustic_tensor(C,QT,theta,N);
        Time_N                =  Time_N + 1;
    end
end
 
  
%--------------------------------------------------------------
% Plot the least of the minors of the acoustic tensor
%--------------------------------------------------------------
switch plotting
    case 1
        q                    =  repmat(q,2,1);
        Nodes_alpha(:,3)     =  q/str.properties.material_parameters.lambda;
        %Nodes_alpha(:,3)    =  q;
        TECPLOT_preprocessing_electro_acoustic_tensor(data,Nodes_alpha,q,connectivity_alpha,1,'acoustic_tensor_simulation_x0');
end



minimum_q                     =  min(q);