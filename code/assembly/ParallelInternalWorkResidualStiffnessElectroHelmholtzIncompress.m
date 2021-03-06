function element_assembly     =  ParallelInternalWorkResidualStiffnessElectroHelmholtzIncompress(ielem,quadrature,solution,geometry,mesh,fem,...
                                                                                           vectorisation,material_information)
                                                                     
%--------------------------------------------------------------------------
% Number of Gauss points
%--------------------------------------------------------------------------
ngauss                        =  size(quadrature.volume.bilinear.Chi,1);
%--------------------------------------------------------------------------
% Initialise assembled residuals per element 
%--------------------------------------------------------------------------
element_assembly              =  ElementResidualMatricesInitialisationElectroIncompressible(geometry,mesh);
%--------------------------------------------------------------------------
% Obtain gradients of kinematics and electrical variables
%--------------------------------------------------------------------------
[kinematics,DN_X_x]           =  KinematicsFunctionVolume(geometry.dim,...
                                          solution.x.Eulerian_x(:,mesh.volume.x.connectivity(:,ielem)),...
                                          solution.x.Lagrangian_X(:,mesh.volume.x.connectivity(:,ielem)),...
                                          fem.volume.bilinear.x.DN_chi);
[E0,DN_X_phi]                 =  ElectricFieldComputation(geometry.dim,...
                                          solution.phi(mesh.volume.phi.connectivity(:,ielem)),...
                                          fem.volume.bilinear.phi.DN_chi,...
                                          kinematics.DX_chi);
%--------------------------------------------------------------------------
% Obtain pressure at every gauss point 
%--------------------------------------------------------------------------
pressure                      =  ScalarFEMInterpolation(fem.volume.bilinear.pressure.N,...
                                               solution.pressure(mesh.volume.pressure.connectivity(:,ielem)));                                         
%--------------------------------------------------------------------------
% First and second derivatives of the model. 
%--------------------------------------------------------------------------
material_information          =  GetDerivativesModelElectroHelmholtz(ielem,geometry.dim,ngauss,kinematics.F,kinematics.H,kinematics.J,E0,...
                                                         material_information,vectorisation);
%--------------------------------------------------------------------------
% First Piola-Kirchhoff stress tensor.        
%--------------------------------------------------------------------------
Piola                         =  FirstPiolaKirchhoffStressTensorElectroHelmholtzIncompressible(ngauss,geometry.dim,kinematics.F,kinematics.H,...
                                                               material_information.derivatives.DPsi.DPsiDF,...
                                                               material_information.derivatives.DPsi.DPsiDH,...
                                                               material_information.derivatives.DPsi.DPsiDJ,...
                                                               pressure);
Piola_vectorised              =  Matrix2Vector(geometry.dim^2,ngauss,Piola);
%--------------------------------------------------------------------------
% Matrix BF    
%--------------------------------------------------------------------------
BF                            =  BMatrix(ngauss,geometry.dim,mesh.volume.x.n_node_elem,DN_X_x,...
                                      vectorisation.Bx_matrix.LHS_indices,...
                                      vectorisation.Bx_matrix.RHS_indices);
%--------------------------------------------------------------------------
% Q matrices arising from the linearisation of H: DH[].SigmaH = Q*DF[]. and
%--------------------------------------------------------------------------
QF                            =  QMatrixComputation(kinematics.F,geometry.dim,ngauss);
QSigmaH                       =  QMatrixComputation(material_information.derivatives.DPsi.DPsiDH,geometry.dim,ngauss);    
%--------------------------------------------------------------------------
% Electric displacement field
%--------------------------------------------------------------------------
D0                            =  -material_information.derivatives.DPsi.DPsiDE0;
%--------------------------------------------------------------------------
% Residuals and Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight        =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.volume.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Residual conservation of linear momentum.
    %----------------------------------------------------------------------
    element_assembly.Tx       =  element_assembly.Tx   +  ...
                                 (BF(:,:,igauss)'*Piola_vectorised(:,igauss))*Integration_weight;
    %----------------------------------------------------------------------
    % Residual Gauss' law.- (D0'*DN_X*W*J_t)'
    %----------------------------------------------------------------------
    element_assembly.Tphi     =  element_assembly.Tphi   -  ...
                                 (DN_X_phi(:,:,igauss)'*D0(:,igauss))*Integration_weight;
    %----------------------------------------------------------------------
    % Residual pressure (incompressibility). 
    %----------------------------------------------------------------------
    element_assembly.Tp       =  element_assembly.Tp  +  ...
                                 fem.volume.bilinear.pressure.N(:,igauss)*(kinematics.J(igauss) - 1)*Integration_weight;
end      
%--------------------------------------------------------------------------
% Residuals and Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight        =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.volume.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices
    %----------------------------------------------------------------------
    vectorised_matrices       =  VectorisedStiffnessMatricesElectroHelmholtzIncompressible(igauss,BF(:,:,igauss),DN_X_phi(:,:,igauss),...
                                                           material_information.derivatives.D2Psi,...
                                                           material_information.derivatives.DPsi,...
                                                           pressure,QF(:,:,igauss),QSigmaH(:,:,igauss),...
                                                           kinematics.H(:,:,igauss),...
                                                           fem.volume.bilinear.pressure.N(:,igauss));
    %----------------------------------------------------------------------
    % Stiffness matrices
    %----------------------------------------------------------------------
    element_assembly.Kxx      =  element_assembly.Kxx      +  vectorised_matrices.Kxx*Integration_weight;
    element_assembly.Kxphi    =  element_assembly.Kxphi    +  vectorised_matrices.Kxphi*Integration_weight;
    element_assembly.Kxp      =  element_assembly.Kxp      +  vectorised_matrices.Kxp*Integration_weight;
    element_assembly.Kphix    =  element_assembly.Kphix    +  vectorised_matrices.Kphix*Integration_weight;
    element_assembly.Kphiphi  =  element_assembly.Kphiphi  +  vectorised_matrices.Kphiphi*Integration_weight;
    element_assembly.Kpx      =  element_assembly.Kpx      +  vectorised_matrices.Kpx*Integration_weight;
    element_assembly.Kpp      =  element_assembly.Kpp      +  vectorised_matrices.Kpp*Integration_weight;
end





