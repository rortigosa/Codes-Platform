%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Residual and Stiffness matrices due to vaccum effect in the principle of
%  virtual work
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function str                 =  BoundaryIntegralsFEM(str)

%--------------------------------------------------------------------------
% Dimension of the problem 
%--------------------------------------------------------------------------
fprintf('Begining of static assembly for boundary integrals in the principle of virtual work\n')
%--------------------------------------------------------------------------
% Dofs per element and initialisation of indexi, indexj and data
%--------------------------------------------------------------------------
% [Kindexi,Kindexj,Kdata,...
% Tindexi,Tindexj,Tdata]      =  SparseStiffnessPreallocationBoundaryFEM(geom,mesh,formulation);
%--------------------------------------------------------------------------
% Loop over elements for the assembly of resiuals and stiffness matrices
%--------------------------------------------------------------------------
tic     
for iedge=1:size(str.mesh.surface.x.boundary_edges,2)
    %----------------------------------------------------------------------
    % Residuals and stiffness matrices
    %----------------------------------------------------------------------
    element_assembly         =  ResidualStiffnessElectroBoundaryFEM(iedge,str.geometry.dim,str.quadrature,str.fem,...
                                               str.mesh,str.solution,str.vectorisation,str.material_information);    
    %----------------------------------------------------------------------
    % Sparse assembly
    %----------------------------------------------------------------------
%     [INDEXI,INDEXJ,DATA]     =  StiffnessSparseAssemblyBoundaryFEM(iedge,str.geometry.dim,str.mesh,element_assembly,str.data.formulation);
%     initial                  =  n_dofs_elem^2*(ielement-1) +  1;
%     final                    =  n_dofs_elem^2*ielement;
%     indexi(initial:final,1)  =  INDEXI;
%     indexj(initial:final,1)  =  INDEXJ;
%     data(initial:final,1)    =  DATA;
%     %----------------------------------------------------------------------
%     % Assembly of residuals
%     %----------------------------------------------------------------------
%     str.assembly             =  VectorsAssemblyBoundaryFEM(iedge,str.geometry.dim,str.mesh,element_assembly,str.assembly,str.data.formulation);
end     
toc 
%--------------------------------------------------------------------------
% Sparse assembly of the stiffness matrix.      
%--------------------------------------------------------------------------
total_dofs                   =  size(str.assembly.K_total,1);
str.assembly.K_total         =  str.assembly.K_total + sparse(indexi,indexj,data,total_dofs,total_dofs);

fprintf('End of static assembly\n')

end





