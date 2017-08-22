function element_assembly  =  ParallelInternalWorkResidualStiffnessUP(ielem,quadrature,solution,geometry,mesh,fem,...
                                                                                           vectorisation,material_information)
ngauss                     =  size(quadrature.volume.bilinear.Chi,1);
%--------------------------------------------------------------------------
% Initialise assembled residuals per element 
%--------------------------------------------------------------------------
element_assembly           =  ElementResidualMatricesInitialisationUP(geometry,mesh);
%--------------------------------------------------------------------------
% Obtain gradients of kinematics and electrical variables
%--------------------------------------------------------------------------
[kinematics,DN_X_x]        =  KinematicsFunctionVolume(geometry.dim,...
                                          solution.x.Eulerian_x(:,mesh.volume.x.connectivity(:,ielem)),...
                                          solution.x.Lagrangian_X(:,mesh.volume.x.connectivity(:,ielem)),...
                                          fem.volume.bilinear.x.DN_chi);
%--------------------------------------------------------------------------
% Obtain pressure at every gauss point 
%--------------------------------------------------------------------------
pressure                   =  ScalarFEMInterpolation(fem.volume.bilinear.pressure.N,...
                                               solution.pressure(mesh.volume.pressure.connectivity(:,ielem)));                                         
%--------------------------------------------------------------------------
% First and second derivatives of the model. 
%--------------------------------------------------------------------------
material_information       =  GetDerivativesModelMechanics(ielem,geometry.dim,ngauss,kinematics.F,kinematics.H,kinematics.J,...
                                                         material_information);
%--------------------------------------------------------------------------
% First Piola-Kirchhoff stress tensor.        
%--------------------------------------------------------------------------
Piola                      =  FirstPiolaKirchhoffStressTensorUP(ngauss,geometry.dim,kinematics.F,kinematics.H,...
                                                               material_information.derivatives.DU.DUDF,...
                                                               material_information.derivatives.DU.DUDH,...
                                                               material_information.derivatives.DU.DUDJ,...
                                                               pressure);
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
    %----------------------------------------------------------------------
    % Residual pressure (incompressibility). 
    %----------------------------------------------------------------------
    element_assembly.Tp    =  element_assembly.Tp  +  ...
                              fem.volume.bilinear.pressure.N(:,igauss)*(kinematics.J(igauss) - 1)*Integration_weight;
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
    vectorised_matrices    =  VectorisedStiffnessMatricesUP(igauss,BF(:,:,igauss),...
                                                      material_information.derivatives.D2U,...
                                                      material_information.derivatives.DU,...
                                                      pressure,QF(:,:,igauss),QSigmaH(:,:,igauss),...
                                                      kinematics.H(:,:,igauss),...
                                                      fem.volume.bilinear.pressure.N(:,igauss));
    %----------------------------------------------------------------------
    % Stiffness matrices
    %----------------------------------------------------------------------
    element_assembly.Kxx   =  element_assembly.Kxx  +  vectorised_matrices.Kxx*Integration_weight;
    element_assembly.Kxp   =  element_assembly.Kxp  +  vectorised_matrices.Kxp*Integration_weight;
    element_assembly.Kpx   =  element_assembly.Kpx  +  vectorised_matrices.Kpx*Integration_weight;
    element_assembly.Kpp   =  element_assembly.Kpp  +  vectorised_matrices.Kpp*Integration_weight;
end




