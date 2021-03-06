%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Sparse assembly for the stiffness matrices arising from the vaccum
%  contribution in the principle of virtual work (boundary integrals)
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [INDEXI,INDEXJ,...
            DATA]           =  StiffnessSparseAssemblyBoundaryBEM(inode,dim,iedge,mesh,...
                                                                element_assembly)   
%--------------------------------------------------------------------------
% x dof's. 
%--------------------------------------------------------------------------
x_nodes                     =  mesh.surface.x.boundary_edges(:,iedge)';
x_dof                       =  zeros(dim,size(x_nodes,2));
for idim=1:dim
    x_dof(idim,:)           =  (x_nodes-1)*dim + idim;
end
x_dof                       =  reshape(x_dof,size(x_dof,1)*size(x_dof,2),1);
n_dofs                      =  mesh.dim*mesh.volume.x.n_nodes; 
%--------------------------------------------------------------------------
% phi dof's. 
%--------------------------------------------------------------------------
phi_dof                     =  mesh.surface.phi.boundary_edges(:,iedge);
phi_dof                     =  n_dofs + phi_dof;
%--------------------------------------------------------------------------
%  phiprime and qprime dof's
%--------------------------------------------------------------------------
if iedge==1
   connectivity             =  mesh.surface.q0.connectivity;
   [local_node,q_elem]      =  find(connectivity==inode);
   q_elem                   =  q_elem(1);
   phiprime_dof             =  mesh.surface.phi.boundary_edges(:,q_elem);
   qprime_dof               =  mesh.surface.q.connectivity(local_node,q_elem);
else
   phiprime_dof             =  []; 
   qprime_dof               =  [];
end    
n_dofs                      =  n_dofs + mesh.volume.phi.n_nodes;
%--------------------------------------------------------------------------
%  pressure dof's
%--------------------------------------------------------------------------
n_dofs                      =  n_dofs + mesh.volume.pressure.n_nodes;
%--------------------------------------------------------------------------
%  q0 dof's
%--------------------------------------------------------------------------
q_dof                       =  mesh.surface.q.connectivity(:,iedge);
q_dof                       =  n_dofs + q_dof;
%--------------------------------------------------------------------------
% Kqx matrix.
%--------------------------------------------------------------------------
newindexi                   =  zeros(size(q_dof,1),size(x_dof,1));
newindexj                   =  zeros(size(q_dof,1),size(x_dof,1));
for iloop=1:size(x_dof,1)
    newindexi(:,iloop)      =  q_dof;
    newindexj(:,iloop)      =  x_dof(iloop)*ones(size(q_dof,1),1);
end
newdata                     =  element_assembly.Kqx(:);
newindexi                   =  reshape(newindexi,size(newindexi,1)*size(newindexi,2),1);
newindexj                   =  reshape(newindexj,size(newindexj,1)*size(newindexj,2),1);
newdata                     =  reshape(newdata,size(newdata,1)*size(newdata,2),1);
Kqx_indexi                  =  newindexi;
Kqx_indexj                  =  newindexj;
Kqx_data                    =  newdata;
%--------------------------------------------------------------------------
% Kqphi matrix.
%--------------------------------------------------------------------------
newindexi                   =  zeros(size(q_dof,1),size(phi_dof,1));
newindexj                   =  zeros(size(q_dof,1),size(phi_dof,1));
for iloop=1:size(phi_dof,1)
    newindexi(:,iloop)      =  q_dof;
    newindexj(:,iloop)      =  phi_dof(iloop)*ones(size(q_dof,1),1);
end
newdata                     =  element_assembly.Kqphi(:);
newindexi                   =  reshape(newindexi,size(newindexi,1)*size(newindexi,2),1);
newindexj                   =  reshape(newindexj,size(newindexj,1)*size(newindexj,2),1);
newdata                     =  reshape(newdata,size(newdata,1)*size(newdata,2),1);
Kqphi_indexi                =  newindexi;
Kqphi_indexj                =  newindexj;
Kqphi_data                  =  newdata;
%--------------------------------------------------------------------------
% Kqq matrix.
%--------------------------------------------------------------------------
newindexi                   =  zeros(size(q_dof,1),size(q_dof,1));
newindexj                   =  zeros(size(q_dof,1),size(q_dof,1));
for iloop=1:size(q_dof,1)
    newindexi(:,iloop)      =  q_dof;
    newindexj(:,iloop)      =  q_dof(iloop)*ones(size(q_dof,1),1);
end
newdata                     =  element_assembly.Kqq(:);
newindexi                   =  reshape(newindexi,size(newindexi,1)*size(newindexi,2),1);
newindexj                   =  reshape(newindexj,size(newindexj,1)*size(newindexj,2),1);
newdata                     =  reshape(newdata,size(newdata,1)*size(newdata,2),1);
Kqq_indexi                  =  newindexi;
Kqq_indexj                  =  newindexj;
Kqq_data                    =  newdata;
%--------------------------------------------------------------------------
% Kqphiprime matrix.
%--------------------------------------------------------------------------
newindexi                   =  zeros(size(qprime_dof,1),size(phiprime_dof,1));
newindexj                   =  zeros(size(qprime_dof,1),size(phiprime_dof,1));
for iloop=1:size(phiprime_dof,1)
    newindexi(:,iloop)      =  q_dof;
    newindexj(:,iloop)      =  phiprime_dof(iloop)*ones(size(q_dof,1),1);
end
newdata                     =  element_assembly.Kqphiprime(:);
newindexi                   =  reshape(newindexi,size(newindexi,1)*size(newindexi,2),1);
newindexj                   =  reshape(newindexj,size(newindexj,1)*size(newindexj,2),1);
newdata                     =  reshape(newdata,size(newdata,1)*size(newdata,2),1);
Kqphiprime_indexi           =  newindexi;
Kqphiprime_indexj           =  newindexj;
Kqphiprime_data             =  newdata;
%--------------------------------------------------------------------------
% Updating indexi, indexj and data. 
%--------------------------------------------------------------------------
INDEXI                      =  [Kqx_indexi;   Kqphi_indexi;   Kqq_indexi;...
                                Kqphiprime_indexi];
INDEXJ                      =  [Kqx_indexj;   Kqphi_indexj;   Kqq_indexj;...
                                Kqphiprime_indexj];
DATA                        =  [Kqx_data;   Kqphi_data;   Kqq_data;...
                                Kqphiprime_data];
end