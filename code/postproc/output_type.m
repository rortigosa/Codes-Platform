function [str]                                     =  output_type(str)


nodes_counter                                      =  str.nodes_counter;
switch str.postproc.output_type
    case 'averaging' 
         %----------------------------------------------------
         % Cauchy stress and different contributions.
         %----------------------------------------------------
         for i = 1:size(str.postproc.stress,1) 
             str.postproc.stress(i,:)              =  str.postproc.stress(i,:)./nodes_counter;
         end

         for i = 1:size(str.postproc.stress,1)
             str.postproc.stress_deviatoric(i,:)   =  str.postproc.stress_deviatoric(i,:)./nodes_counter;
         end

         for i = 1:size(str.postproc.stress,1)
             str.postproc.stress_volumetric(i,:)   =  str.postproc.stress_volumetric(i,:)./nodes_counter;
         end
         %----------------------------------------------------
         % Stress conjugates.
         %----------------------------------------------------
         for i = 1:size(str.postproc.SigmaF,1) 
             for j = 1:size(str.postproc.SigmaF,2) 
                 str.postproc.SigmaF(i,j,:)        =  reshape(str.postproc.SigmaF(i,j,:),size(str.postproc.SigmaF(i,j,:),3),1)'./nodes_counter;
             end
         end

         for i = 1:size(str.postproc.SigmaH,1) 
             for j = 1:size(str.postproc.SigmaH,2) 
             str.postproc.SigmaH(i,j,:)            =  reshape(str.postproc.SigmaH(i,j,:),size(str.postproc.SigmaH(i,j,:),3),1)'./nodes_counter;
             end
         end

         for i = 1:str.postproc.n_nodes
             str.postproc.SigmaJ(i,1)              =  str.postproc.SigmaJ(i,1)/nodes_counter(i);
         end  
         
         for i = 1:size(str.postproc.SigmaF,1) 
             for j = 1:size(str.postproc.SigmaF,2) 
             str.postproc.First_Piola(i,j,:)       =  reshape(str.postproc.First_Piola(i,j,:),size(str.postproc.First_Piola(i,j,:),3),1)'./nodes_counter;
             end
         end         
         %----------------------------------------------------
         % Cauchy and different contributions.
         %----------------------------------------------------
         for i = 1:size(str.postproc.stress,1)
             str.postproc.sigma(i,:)               =  str.postproc.sigma(i,:)./nodes_counter;
         end

         for i = 1:size(str.postproc.sigma,1)
             str.postproc.sigma_deviatoric(i,:)    =  str.postproc.sigma_deviatoric(i,:)./nodes_counter;
         end

         for i = 1:size(str.postproc.sigma,1)
             str.postproc.sigma_volumetric(i,:)    =  str.postproc.sigma_volumetric(i,:)./nodes_counter;
         end

         str.postproc.sigma_pressure               =  str.postproc.sigma_pressure'./nodes_counter;
         %----------------------------------------------------
         % Lagrangian Electric displacement.
         %----------------------------------------------------
         for i = 1:size(str.postproc.D0,1)
             str.postproc.D0(i,:)                  =  str.postproc.D0(i,:)./nodes_counter;
         end
        
         %----------------------------------------------------
         % Eulerian Electric displacement.
         %----------------------------------------------------
         for i = 1:size(str.postproc.D,1)
             str.postproc.D(i,:)                   =  str.postproc.D(i,:)./nodes_counter;
         end
         %----------------------------------------------------
         % Electric field.
         %----------------------------------------------------
         for i = 1:size(str.postproc.E0,1)
             str.postproc.E0(i,:)                  =  str.postproc.E0(i,:)./nodes_counter;
         end

         for i = 1:size(str.postproc.E,1)
             str.postproc.E(i,:)                   =  str.postproc.E(i,:)./nodes_counter;
         end
         str.postproc.E_norm                       =  str.postproc.E_norm'./nodes_counter;
         str.postproc.electric_breakdown_factor    =  str.postproc.electric_breakdown_factor'./nodes_counter;
         %----------------------------------------------------
         % Deformation variables.
         %----------------------------------------------------
         for i = 1:str.postproc.n_nodes
             str.postproc.F(:,:,i)                 =  str.postproc.F(:,:,i)/nodes_counter(i);
         end
         for i = 1:str.postproc.n_nodes
             str.postproc.H(:,:,i)                 =  str.postproc.H(:,:,i)/nodes_counter(i);
         end
         for i = 1:str.postproc.n_nodes
             str.postproc.J(i,1)                   =  str.postproc.J(i,1)/nodes_counter(i);
         end         
         for i = 1:str.postproc.n_nodes
             str.postproc.FL2norm(i,1)             =  str.postproc.FL2norm(i,1)/nodes_counter(i);
         end         
         switch str.data.formulation.mixed_type
               case {'twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v3_mixed_formulation','twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
                    for i = 1:size(str.postproc.E,1)
                        str.postproc.gradJ(i,:)        =  str.postproc.gradJ(i,:)./nodes_counter;
                    end
                    for i = 1:str.postproc.n_nodes
                        str.postproc.ratioJ(i,1)        =  str.postproc.ratioJ(i,1)/nodes_counter(i);
                    end
         end                   
         %----------------------------------------------------
         % Plastic deformation.
         %----------------------------------------------------
         for i = 1:str.postproc.n_nodes
             str.postproc.hardening_variable(i)    =  str.postproc.hardening_variable(i)/nodes_counter(i);
             str.postproc.lambda_p(:,i)            =  str.postproc.lambda_p(:,i)/nodes_counter(i);
         end
         %----------------------------------------------------
         % Intensity density vector.
         %----------------------------------------------------
         for i = 1:size(str.postproc.J_electric,1)
             str.postproc.J_electric(i,:)          =  str.postproc.J_electric(i,:)./nodes_counter;
         end
         
         for i = 1:size(str.postproc.J0_electric,1)
             str.postproc.J0_electric(i,:)         =  str.postproc.J0_electric(i,:)./nodes_counter;
         end
         %----------------------------------------------------
         % Eigenvalues of the second derivative of the internal
         % energy with respect to the deformation gradient F.
         %----------------------------------------------------
         for i = 1:size(str.postproc.eigenvalues,1)
             str.postproc.eigenvalues(i,:)                  =  str.postproc.eigenvalues(i,:)./nodes_counter;
         end

end