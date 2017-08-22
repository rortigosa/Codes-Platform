%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function computes the elemental residual vectors and stiffness
% matrices for the formulation with fields: x-phi-p
% 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                =  ResidualStiffnessElectroBoundaryBEM(inode,iedge,str)

%--------------------------------------------------------------------------
% Coordinate of the collocation point and the Gauss points in the boundary
%--------------------------------------------------------------------------
switch str.mesh.dim
    case 2
        dim_factor          =  2*pi;
    case 3
        dim_factor          =  -4*pi;
end
%--------------------------------------------------------------------------
% Coordinate of the collocation point and the Gauss points in the boundary
%--------------------------------------------------------------------------
xedge                       =  str.mesh.surface.x.boundary_edges(:,iedge);
xGauss                      =  vector_FEM_interpolation(str.fem.surface.BEM_FEM.x.N,xedge);
%--------------------------------------------------------------------------
% Find one of the elements that xprime belongs to and its local numbering
% within that element
%--------------------------------------------------------------------------
connectivity               =  str.mesh.surface.q0.connectivity;
[local_node,q0_elem]       =  find(connectivity==inode);
local_node                 =  local_node(1);
q0_elem                    =  q0_elem(1);
xprime_elem                =  str.solution.x.Eulerian_x(:,str.mesh.surface.x.boundary_edges(:,q0_elem));
xprime                     =  vector_FEM_interpolation(str.fem.surface.nodes.x.N_q0(:,local_node),xprime_elem);
%--------------------------------------------------------------------------
% Distance between xGauss and xprime
%--------------------------------------------------------------------------
r                          =  xGauss - repmat(xprime,size(xGauss,2));
r_norm                     =  VectorNorm(r); 
%--------------------------------------------------------------------------
% Test functions in BEM 
%--------------------------------------------------------------------------
V                           =  Laplace_Fundamental_Solution(str.mesh.dim,r_norm);
dVdx                        =  diff_Laplace_Fundamental_Solution(str.mesh.dim,r,r_norm);
d2Vdxdx                     =  diff_diff_Laplace_Fundamental_Solution(str.mesh.dim,r,r_norm);
%--------------------------------------------------------------------------
% Electric potential 
%--------------------------------------------------------------------------
phiedge                     =  str.mesh.surface.phi.boundary_edges(:,iedge);
phiGauss                    =  scalar_FEM_interpolation(str.fem.surface.BEM_FEM.phi.N,phiedge);
phiprime_elem               =  str.solution.x.Eulerian_x(:,str.mesh.surface.x.boundary_edges(:,q0_elem));
phiprime                    =  vector_FEM_interpolation(str.fem.surface.nodes.phi.N_q0(:,local_node),phiprime_elem);
%--------------------------------------------------------------------------
% Obtain the value of the flux field q0
%--------------------------------------------------------------------------
q0                          =  scalar_FEM_interpolation(fem.surface.bilinear.q0.N,solution.q0(sr.mesh.surface.q0.connectivity(:,iedge)));
%--------------------------------------------------------------------------
% Obtain gradients of kinematics and electrical variables
%--------------------------------------------------------------------------
[kinematics,DN_X_x]         =  kinematics_function_surface(str.mesh.dim,...
                                         str.solution.x.Eulerian_x(:,str.mesh.surface.x.boundary_edges(:,ielem)),...
                                         str.solution.x.Lagrangian_X(:,str.mesh.surface.x.boundary_edges(:,ielem)),...
                                         str.solution.x.Eulerian_x(:,str.mesh.volume.x.connectivity(:,ielem)),...
                                         str.fem.surface.BEM_FEM.x.DN_chi);
                                         
H_vectorisedT               =  reshape(permute(kinematics.H,[2 1 3]),1,str.mesh.dim^2,[]);                                         
HN                          =  MatrixVectorMultiplication(kinematics.H,kinematics.Normal_vector);
HN_norm                     =  VectorNorm(HN);
%--------------------------------------------------------------------------
% Matrix BF   
%--------------------------------------------------------------------------
BF                          =  Bmatrix(str.mesh.dim,str.mesh.volume.x.n_node_elem,DN_X_x,...
                                        str.vectorisation.Bx_matrix_boundary_BEM.LHS_indices,...
                                        str.vectorisation.Bx_matrix_boundary_BEM.RHS_indices);
%--------------------------------------------------------------------------
% Required vectorisations 
%--------------------------------------------------------------------------
N_matrix                    =  vector2matrix_vectorisation(kinematics.Normal_vector,n_gauss,str.vectorisation.vector2matrix_rowwise_Gauss_points_boundary_BEM);
N_matrixT                   =  permute(N_matrix,[2 1 3]);
QF                          =  Q_matrix_computation(kinematics.F);
Nx_vectorised               =  str.fem.surface.bilinear.x.N_vectorised;
%--------------------------------------------------------------------------
% Residuals 
%--------------------------------------------------------------------------
for igauss=1:size(str.quadrature.surface.bilinear.Chi,1)
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight      =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.surface.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Residual for q.
    %----------------------------------------------------------------------
    vector1                 =  phiGauss(igauss)*(dVdx(:,igauss)'*HN(:,igauss));
    vector2                 =  -V(igauss)*q0(igauss)*HN_norm(igauss);
    element_assembly.Tq     =  element_assembly.Tq   +  (vector1 + vector2)*Integration_weight;
end 
%--------------------------------------------------------------------------
% Residual Tqprime
%--------------------------------------------------------------------------
if iedge==1
   str.assembly.element_assembly.Tqprime     =  -dim_factor*phiprime;
end
     
%--------------------------------------------------------------------------
% Stiffness matrices 
%--------------------------------------------------------------------------
for igauss=1:size(str.quadrature.volume.bilinear.Chi,1)
    %----------------------------------------------------------------------
    % Integration weight
    %----------------------------------------------------------------------
    Integration_weight      =  (kinematics.DX_chi_Jacobian(igauss))*quadrature.surface.bilinear.W_v(igauss);
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kqx
    %----------------------------------------------------------------------
    vector                  =  N_matrix(:,:,igauss)*QF(:,:,igauss)*BF(:,:,igauss);
    Kqx1                    =  phiGauss(igauss)*HN(:,igauss)'*d2Vdxdx(:,:,igauss)*Nx_vectorised(:,:,igauss);
    Kqx2                    =  phiGauss(igauss)*dVdx(:,igauss)'*vector;
    Kqx3                    =  -q0(igauss)*HN_norm(igauss)*(dVdx(:,igauss)'*Nx_vectorised(:,:,igauss));
    Kqx4                    =  -V(igauss)*q0(igauss)/HN_norm(igauss)*(H_vectorisedT(:,igauss)*N_matrixT(:,:,igauss)*vector);
    Kqx                     =  Kqx1 + Kqx2 + Kqx3 + Kqx4;
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kqphi
    %----------------------------------------------------------------------
    Kqphi                   =  (dVdx(:,igauss)'*HN(:,igauss))*str.fem.surface.BEM_FEM.phi.N(:,igauss);
    %----------------------------------------------------------------------
    % Vectorisation of stiffness matrices Kqq
    %----------------------------------------------------------------------
    Kqq                     =  -V(igauss)*HN_norm(igauss)*str.fem.surface.BEM_FEM.q0.N(:,igauss);
    %----------------------------------------------------------------------
    % Stiffness matrices
    %----------------------------------------------------------------------
    str.assembly.element_assembly.Kqx    =  str.assembly.element_assembly.Kqq     +  Kqx*Integration_weight;
    str.assembly.element_assembly.Kqphi  =  str.assembly.element_assembly.Kqphi   +  Kqphi*Integration_weight;
    str.assembly.element_assembly.Kqq    =  str.assembly.element_assembly.Kqq     +  Kqq*Integration_weight;    
end
%--------------------------------------------------------------------------
% Stiffness matrix Kqphiprime
%--------------------------------------------------------------------------
if iedge==1
   str.assembly.element_assembly.Kqphiprime  =  -dim_factor*str.fem.surface.nodes.phi.N_q0(:,local_node);
end
