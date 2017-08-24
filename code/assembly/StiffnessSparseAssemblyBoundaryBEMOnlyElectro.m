%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Sparse assembly for the stiffness matrices arising from the vaccum
%  contribution in the principle of virtual work (boundary integrals)
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [INDEXI,INDEXJ,...
            DATA]           =  StiffnessSparseAssemblyBoundaryBEMOnlyElectro(inode,dim,iedge,mesh,...
                                                                element_assembly)   



%--------------------------------------------------------------------------
% phi dof's. 
%--------------------------------------------------------------------------
phi_dof                     =  mesh.surface.phi.boundary_edges(:,iedge);
%--------------------------------------------------------------------------
%  phiprime dof's
%--------------------------------------------------------------------------
% if iedge==1
%    %connectivity            =  mesh.surface.q.connectivity(:,iedge);
%    %[local_node,q_elem]     =  find(connectivity==inode);

%    %q_elem                  =  q_elem(1);
%    %phiprime_dof            =  mesh.surface.phi.boundary_edges(:,q_elem);   
%    %qprime_dof              =  mesh.surface.q.connectivity(local_node,q_elem);   
% connectivity        =  mesh.surface.q.connectivity;
% q_elem              =  ceil(find(connectivity'==inode)/dim);
% local_node          =  find(mesh.surface.q.connectivity(:,q_elem(1))==inode);
%    phiprime_dof             =  mesh.surface.phi.boundary_edges(:,q_elem(1));
%    qprime_dof               =  mesh.surface.q.connectivity(local_node,q_elem(1));
% else
%    phiprime_dof             =  []; 
%    qprime_dof               =  []; 
% end    
n_dofs                      =  mesh.volume.phi.n_nodes;
%--------------------------------------------------------------------------
%  q dof's
%--------------------------------------------------------------------------
q_dof                       =  mesh.surface.q.connectivity(:,iedge);
q_dof                       =  n_dofs + q_dof;
qprime_dof                  =  inode + n_dofs;
%--------------------------------------------------------------------------
% Kqphi matrix.
%--------------------------------------------------------------------------
newindexi                   =  zeros(size(qprime_dof,1),size(phi_dof,1));
newindexj                   =  zeros(size(qprime_dof,1),size(phi_dof,1));
for iloop=1:size(phi_dof,1)
    newindexi(:,iloop)      =  qprime_dof;
    newindexj(:,iloop)      =  phi_dof(iloop)*ones(size(qprime_dof,1),1);
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
newindexi                   =  zeros(size(qprime_dof,1),size(q_dof,1));
newindexj                   =  zeros(size(qprime_dof,1),size(q_dof,1));
for iloop=1:size(q_dof,1)
    newindexi(:,iloop)      =  qprime_dof;
    newindexj(:,iloop)      =  q_dof(iloop)*ones(size(qprime_dof,1),1);
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
% newindexi                   =  zeros(size(qprime_dof,1),size(phiprime_dof,1));
% newindexj                   =  zeros(size(qprime_dof,1),size(phiprime_dof,1));
% for iloop=1:size(phiprime_dof,1)
%     newindexi(:,iloop)      =  qprime_dof;
%     newindexj(:,iloop)      =  phiprime_dof(iloop)*ones(size(qprime_dof,1),1);
% end
% newdata                     =  element_assembly.Kqphiprime(:);
% newindexi                   =  reshape(newindexi,size(newindexi,1)*size(newindexi,2),1);
% newindexj                   =  reshape(newindexj,size(newindexj,1)*size(newindexj,2),1);
% newdata                     =  reshape(newdata,size(newdata,1)*size(newdata,2),1);
% Kqphiprime_indexi           =  newindexi;
% Kqphiprime_indexj           =  newindexj;
% Kqphiprime_data             =  newdata;
%--------------------------------------------------------------------------
% Updating indexi, indexj and data. 
%--------------------------------------------------------------------------
INDEXI                      =  [Kqphi_indexi;   Kqq_indexi];
INDEXJ                      =  [Kqphi_indexj;   Kqq_indexj];
DATA                        =  [Kqphi_data;     Kqq_data];
% INDEXI                      =  [Kqphi_indexi;   Kqq_indexi;  Kqphiprime_indexi];
% INDEXJ                      =  [Kqphi_indexj;   Kqq_indexj;  Kqphiprime_indexj];
% DATA                        =  [Kqphi_data;     Kqq_data;    Kqphiprime_data];
end

