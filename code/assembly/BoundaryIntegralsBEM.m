%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Residual and Stiffness matrices due to vaccum effect in the principle of
%  virtual work
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function str                 =  BoundaryIntegralsBEM(str)

%--------------------------------------------------------------------------
% Dimension of the problem 
%--------------------------------------------------------------------------
fprintf('Begining of static assembly for BEM integral\n')
%--------------------------------------------------------------------------
% Initialise assembled residuals per element
%--------------------------------------------------------------------------
str                          =  ElementResidualInitialisationFormulationBoundaryBEM(str);
%--------------------------------------------------------------------------
% Dofs per element and initialisation of indexi, indexj and data
%--------------------------------------------------------------------------
[indexi,indexj,data]         =  SparseStiffnessPreallocationBoundaryBEM(str.mesh);
%--------------------------------------------------------------------------
% Loop over collocation points
%--------------------------------------------------------------------------
tic     
for inode=1:str.mesh.surface.q0.n_nodes
    for iedge=1:size(str.mesh.surface.x.boundary_edges,2)
        %------------------------------------------------------------------
        % Residuals and stiffness matrices
        %------------------------------------------------------------------
        element_assembly         =  ResidualStiffnessElectroBoundaryBEM(inode,iedge,str);
        %------------------------------------------------------------------
        % Sparse assembly
        %------------------------------------------------------------------
        [INDEXI,INDEXJ,DATA]     =  StiffnessSparseAssemblyBoundaryBEM(iedge,str.geometry.dim,str.mesh,element_assembly,str.data.formulation);
        initial                  =  n_dofs_elem^2*(ielement-1) +  1;
        final                    =  n_dofs_elem^2*ielement;
        indexi(initial:final,1)  =  INDEXI;
        indexj(initial:final,1)  =  INDEXJ;
        data(initial:final,1)    =  DATA;
        %------------------------------------------------------------------
        % Assembly of residuals
        %------------------------------------------------------------------
        str.assembly             =  VectorsAssemblyBoundaryBEM(iedge,str.mesh,element_assembly,str.assembly);
    end
end     
toc 
%--------------------------------------------------------------------------
% Sparse assembly of the stiffness matrix.      
%--------------------------------------------------------------------------
total_dofs                   =  size(str.assembly.K_total,1);
str.assembly.K_total         =  str.assembly.K_total + sparse(indexi,indexj,data,total_dofs,total_dofs);

fprintf('End of static assembly\n')

end





