function [INDEXI,...
  INDEXJ,DATA]   =  ParallelInternalVectorsAssemblyOnlyElectro(iedge,mesh,element_assembly)

global_phi_dof   =  mesh.surface.phi.boundary_edges(:,iedge);
global_q_dof     =  mesh.surface.q.connectivity(:,iedge);
%--------------------------------------------------------------------------
% INDEXI, INDEXJ, DATA
%--------------------------------------------------------------------------
INDEXI           =  [global_phi_dof;global_q_dof];
INDEXJ           =  ones(size(INDEXI,1),1);
DATA             =  [element_assembly.Tphi;element_assembly.Tq];
                                    