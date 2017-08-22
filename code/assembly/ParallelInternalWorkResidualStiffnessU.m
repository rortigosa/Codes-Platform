function element_assembly  =  ParallelInternalWorkResidualStiffnessU(ielem,quadrature,solution,geometry,mesh,fem,...
                                                                                           vectorisation,material_information)
ngauss                     =  size(quadrature.volume.bilinear.Chi,1);
%--------------------------------------------------------------------------
% Initialise assembled residuals per element 
%--------------------------------------------------------------------------
element_assembly           =  ElementResidualMatricesInitialisationU(geometry,mesh);
%--------------------------------------------------------------------------
% Obtain gradients of kinematics and electrical variables
%--------------------------------------------------------------------------
[kinematics,DN_X_x]        =  KinematicsFunctionVolume(geometry.dim,...
                                          solution.x.Eulerian_x(:,mesh.volume.x.connectivity(:,ielem)),...
                                          solution.x.Lagrangian_X(:,mesh.volume.x.connectivity(:,ielem)),...
                                          fem.volume.bilinear.x.DN_chi);
%--------------------------------------------------------------------------
% First and second derivatives of the model. 
%--------------------------------------------------------------------------
material_information       =  GetDerivativesModelMechanics(ielem,geometry.dim,ngauss,kinematics.F,kinematics.H,kinematics.J,...
                                                         material_information);
%--------------------------------------------------------------------------
% First Piola-Kirchhoff stress tensor.        
%--------------------------------------------------------------------------
Piola                      =  FirstPiolaKirchhoffStressTensorU(ngauss,geometry.dim,kinematics.F,kinematics.H,...
                                                               material_information.derivatives.DU.DUDF,...
                                                               material_information.derivatives.DU.DUDH,...
                                                               material_information.derivatives.DU.DUDJ);
Piola_vectorised           =  Matrix2Vector(geometry.dim^2,ngauss,Piola);
%--------------------------------------------------------------------------
% Matrix BF    
%--------------------------------------------------------------------------
BF                         =  BMatrix(ngauss,geometry.dim,mesh.volume.x.n_node_elem,DN_X_x,...
                                      vectorisation.Bx_matrix.LHS_indices,...
                                      vectorisation.Bx_matrix.RHS_indices);
%--------------------------------------------------------------------------
% Q matrices arising from the linearisation of H: DH[].SigmaH = Q*DF[]. and
%--------------------------------------------------------------------------
QF                         =  QMatrixComputation(kinematics.F,geometry.dim,ngauss);
QSigmaH                    =  QMatrixComputation(material_information.derivatives.DU.DUDH,geometry.dim,ngauss);    
%--------------------------------------------------------------------------
% Residuals and Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight     =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.volume.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Residual conservation of linear momentum.
    %----------------------------------------------------------------------
    element_assembly.Tx    =  element_assembly.Tx   +  ...
                              (BF(:,:,igauss)'*Piola_vectorised(:,igauss))*Integration_weight;
end      
%--------------------------------------------------------------------------
% Residuals and Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight     =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.volume.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices
    %----------------------------------------------------------------------
    vectorised_matrices    =  VectorisedStiffnessMatricesU(igauss,BF(:,:,igauss),...
                                                      material_information.derivatives.D2U,...
                                                      material_information.derivatives.DU,...
                                                      QF(:,:,igauss),QSigmaH(:,:,igauss),...
                                                      kinematics.H(:,:,igauss));
    %----------------------------------------------------------------------
    % Stiffness matrices
    %----------------------------------------------------------------------
    element_assembly.Kxx   =  element_assembly.Kxx  +  vectorised_matrices.Kxx*Integration_weight;
end




