%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function computes the elemental residual vectors and stiffness
% matrices for the formulation with fields: x-phi-p
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function asmb       =  ResidualStiffnessElectroBoundaryFEM(iedge,dim,quadrature,fem,mesh,solution,...
                                                                   vectorisation,mat_info)
%--------------------------------------------------------------------------
% Number of Gauss points
%--------------------------------------------------------------------------
ngauss              =  size(quadrature.surface.bilinear.Chi,1);
%--------------------------------------------------------------------------
% Permittivity of vaccum
%--------------------------------------------------------------------------
e0                  =  mat_info.vacuum_electric_permittivity;
%--------------------------------------------------------------------------
% Obtain the value of the flux field q0
%--------------------------------------------------------------------------
q                   =  ScalarFEMInterpolation(fem.surface.bilinear.q.N,solution.q(mesh.surface.q.connectivity(:,iedge)));
%---------------------------------------------------------- ----------------
% Obtain gradients of kinematics and electrical variables
%--------------------------------------------------------------------------
[kinematics,DN_X_x] =  KinematicsFunctionSurface(dim,solution.x.Lagrangian_X(:,mesh.surface.x.boundary_edges(:,iedge)),...
                                solution.x.Lagrangian_X(:,mesh.surface.x.boundary_edges(:,iedge)),...
                                solution.x.Lagrangian_X(:,mesh.volume.x.connectivity(:,mesh.surface.x.volume_elements(iedge))),...
                                fem.surface.bilinear.x.N,fem.surface.bilinear.x.DN_chi);                                         
[E0,DN_X_phi]       =  ElectricFieldSurfaceComputation(dim,solution.phi(mesh.surface.x.boundary_edges(:,iedge)),...
                                        fem.surface.bilinear.phi.DN_chi,kinematics.DX_chi);
H_vectorisedT       =  reshape(permute(kinematics.H,[2 1 3]),1,dim^2,[]);                                         
HN                  =  MatrixVectorMultiplication(dim,ngauss,kinematics.H,kinematics.Normal_vector);
HN_norm             =  VectorNorm(HN);
%--------------------------------------------------------------------------
% Spatial electric field
%--------------------------------------------------------------------------
E                   =  MatrixVectorMultiplication(dim,ngauss,kinematics.FminusT,E0);
%--------------------------------------------------------------------------
% Matrix BF   
%--------------------------------------------------------------------------
BF                  =  BMatrix(ngauss,dim,size(mesh.surface.x.boundary_edges,1),DN_X_x,...
                               vectorisation.Bx_matrix_boundary.LHS_indices,...
                               vectorisation.Bx_matrix_boundary.RHS_indices);
%--------------------------------------------------------------------------
% Required vectorisations 
%--------------------------------------------------------------------------
N_matrix            =  Vector2MatrixVectorisation(kinematics.Normal_vector,ngauss,vectorisation.vector2matrix_rowwise_Gauss_points_boundary_FEM);
N_matrixT           =  permute(N_matrix,[2 1 3]);
Finv                =  permute(kinematics.FminusT,[2 1 3]);
invCE0              =  MatrixVectorMultiplication(Finv,E);
invCE0_matrix       =  vector2matrix_vectorisation(invCE0,n_gauss,vectorisation.vector2matrix_rowwise_Gauss_points_boundary_FEM);
B                   =  BmatrixBoundary(E,Finv,vectorisation.DFminusT.LHS_indices,vectorisation.DFminusT.RHS_indices);
BT                  =  permute(B,[2 1 3]);
QF                  =  Q_matrix_computation(kinematics.F);
Nx_vectorised       =  fem.surface.bilinear.x.N_vectorised;
Nx_vectorisedT      =  permute(Nx_vectorised,[2 1 3]);
%--------------------------------------------------------------------------
% Residuals 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Int_weight      =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.surface.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Residual conservation of linear momentum.
    %----------------------------------------------------------------------
    Tx1             =  -e0*q0(igauss)*norm(HN_norm(igauss))*E(:,igauss);
    Tx2             =  e0/2*(E(:,igauss)'*E(:,igauss) - q0(igauss)^2)*HN(:,igauss);    
    asmb.Tx         =  asmb.Tx   +  Nx_vectorisedT(:,:,igauss)*(Tx1 + Tx2)*Int_weight;
    %----------------------------------------------------------------------
    % Residual Gauss' law.- (D0'*DN_X*W*J_t)'
    %----------------------------------------------------------------------
    asmb.Tphi       =  asmb.Tphi   -  (q0*HN_norm(igauss)*fem.surface.bilinear.phi.N(:,igauss))*Int_weight;
end      
%--------------------------------------------------------------------------
% Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:ngauss
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Int_weight      =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.surface.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kxx
    %----------------------------------------------------------------------
    vector1         =  Nx_vectorisedT*E(:,igauss);
    vector2         =  (H_vectorisedT(:,igauss)*N_matrixT(:,:,igauss))*(N_matrix(:,:,igauss)*QF(:,:,igauss)*BF(:,:,igauss));
    Kxx1            =  -e0*q0(igauss)/HN_norm(igauss)*(vector1*vector2');
    
    Kxx2            =  e0*q0(igauss)*HN_norm(igauss)*(N_matrixT(:,:,igauss)*BT(:,:,igauss)*BF(:,:,igauss));
    
    Kxx3            =  -e0/2*q0(igauss)^2*(N_matrixT(:,:,igauss)*N_matrix(:,:,igauss)*QF(:,:,igauss));
    
    vector3         =  Nx_vectorisedT(:,:,igauss)*HN(:,igauss);
    vector4         =  E(:,igauss)'*invCE0_matrix(:,:,igauss)*BF(:,:,igauss);
    Kxx4            =  -e0*(vector3*vector4');
    
    Kxx5            =  e0*(E(:,igauss)'*E(:,igauss))*(Nx_vectorisedT(:,:,igauss)*N_matrix(:,:,igauss)*QF(:,:,igauss)*BF(:,:,igauss));
    
    Kxx             =  Kxx1 + Kxx2 + Kxx3 + Kxx4 + Kxx5;
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kxphi
    %----------------------------------------------------------------------
    vector5         =  Nx_vectorisedT(:,:,igauss)*kinematics.FminusT(:,:,igauss)*DN_X_phi(:,:,igauss);
    Kxphi1          =  -e0*q0(igauss)*HN_norm(igauss)*vector5;
    
    vector6         =  (Nx_vectorisedT(:,:,igauss)*HN(:,igauss));
    vector7         =  invCE0(:,igauss)'*DN_X_phi(:,:,igauss);
    Kxphi2          =  e0*q0(igauss)*(vector6*vector7);
    
    Kxphi           =  Kxphi1 + Kxphi2;
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kxq
    %----------------------------------------------------------------------
    Kxq1            =  -e0*HN_norm(igauss)*(N_matrixT(:,:,igauss)*E(:,igauss))*fem.surface.bilinear.q0.N(:,igauss)';
    Kxq2            =  e0*q0(igauss)*(N_matrixT(:,:,igauss)*HN(:,igauss))*fem.surface.BEM_FEM.q0.N(:,igauss)';
    Kxq             =  Kxq1 + Kxq2;
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kphix
    %----------------------------------------------------------------------
    vector8         =  H_vectorisedT(:,igauss)*N_matrixT(:,:,igauss)*N_matrix(:,:,igauss)*QF(:,:,igauss)*BF(:,:,igauss);
    Kphix           =  -(e0*q0(igauss)/HN_norm(igauss))*fem.surface.bilinear.phi.N(:,igauss)*vector8;
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kphiq
    %----------------------------------------------------------------------
    Kphiq           =  -e0*HN_norm(igauss)*(fem.surface.bilinear.phi.N(:,gauss)*fem.surface.BEM_FEM.q0.N(:,igauss)');
    %----------------------------------------------------------------------
    % Stiffness matrices
    %----------------------------------------------------------------------
    asmb.Kxx      =  asmb.Kxx      +  Kxx*Int_weight;
    asmb.Kxphi    =  asmb.Kxphi    +  Kxphi*Int_weight;
    asmb.Kxq      =  asmb.Kxq      +  Kxq*Int_weight;    
    asmb.Kphix    =  asmb.Kphix    +  Kphix*Int_weight;
    asmb.Kphip    =  asmb.Kphip    +  Kphiq*Int_weight;
end