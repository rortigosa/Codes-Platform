%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Assemblign of the resiual and stiffness matrix of the system.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                 =  MassMatricesAssembly(str)    

%--------------------------------------------------------------------------
% Dimension of the problem 
%--------------------------------------------------------------------------
fprintf('Begining of mass matrices assembly\n')
%--------------------------------------------------------------------------
% Dofs per element and initialisation of indexi, indexj and data
%--------------------------------------------------------------------------
n_dofs_elem_x0               =  str.mesh.volume.x.n_node_elem*3;

n_dofs_elem                  =  n_dofs_elem_x0^2;
n_dofs                       =  n_dofs_elem*str.mesh.volume.n_elem;
indexi                       =  zeros(n_dofs,1);
indexj                       =  zeros(n_dofs,1);
data                         =  zeros(n_dofs,1);
%--------------------------------------------------------------------------
% Loop over elements for the assembly of resiuals and stiffness matrices
%--------------------------------------------------------------------------
for ielement=1:str.mesh.volume.n_elem 
    %----------------------------------------------------------------------
    % Sparse assembly
    %----------------------------------------------------------------------
    element_assembly         =  MassMatricesComputation(ielement,str.geometry.dim,str.solution,str.mesh,str.quadrature.volume.mass,...
                                                          str.fem.volume.mass,str.assembly.mass_matrices,str.material_information);    
    [INDEXI,INDEXJ,DATA]     =  MassSparseAssembly(ielement,str.geometry.dim,str.mesh,element_assembly);
    
    initial                  =  n_dofs_elem*(ielement-1) + 1;
    final                    =  n_dofs_elem*(ielement);
    indexi(initial:final,1)  =  INDEXI;
    indexj(initial:final,1)  =  INDEXJ;
    data(initial:final,1)    =  DATA;
end   
%--------------------------------------------------------------------------
% Sparse assembly of the stiffness matrix.      
%--------------------------------------------------------------------------
total_dofs                   =  str.solution.n_dofs;
str.assembly.M_total         =  sparse(indexi,indexj,data,total_dofs,total_dofs);

fprintf('End of mass matrices assembly\n')

end





