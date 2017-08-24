%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Residual arising from the BEM integral
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [INDEXI,INDEXJ,...
    DATA]                  =  SparseVectorsAssemblyBoundaryBEM(iedge,mesh,element_assembly,assembly)

global_q_dof               =  mesh.surface.q.connectivity(:,iedge) + ;
DATA                       =  element_assembly.Tq(:);         

