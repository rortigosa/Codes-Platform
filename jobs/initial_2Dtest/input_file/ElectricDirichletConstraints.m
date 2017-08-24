%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  In this function we specify boundary conditions for the mechanical
%  physics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    [fixdof_phi,freedof_phi,...
    cons_val_phi]        =  ElectricDirichletConstraints(str)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Dirichlet boundary conditions for displacements (I will specify that all 
% the degrees of freedom at z=0 are fixed. A prescribed displacement at 
% z=1 of Dz=0.1 will be imposed in the Z direction)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Displacements 
%--------------------------------------------------------------------------
nodes1                   =  zeros(str.mesh.volume.phi.n_nodes,1);
presc1                   =  zeros(str.mesh.volume.phi.n_nodes,1);
ymed                     =  (min(str.mesh.volume.phi.nodes(2,:)) + max(str.mesh.volume.phi.nodes(2,:)))/2;
for inode=1:str.mesh.volume.phi.n_nodes
    x                    =  str.mesh.volume.phi.nodes(:,inode);
    if (x(2) - ymed)<1e-6
       nodes1(inode)     =  inode;
       presc1(inode)     =  0;
    end 
end
nodes1                   =  nodes1(nodes1>0);
presc1                   =  presc1(nodes1);
%--------------------------------------------------------------------------
% Degrees of freedom with associated constraints
%--------------------------------------------------------------------------
dofs1                    =  nodes1;  %  Select direction X for set 1
freedof_phi              =  (1:str.mesh.volume.phi.n_nodes)';
[fixdof_phi,order]       =  sort(dofs1);
freedof_phi(fixdof_phi)  =  [];
%--------------------------------------------------------------------------
% Degrees of freedom with associated constraints 
%--------------------------------------------------------------------------
prescribed               =  presc1;
prescribed               =  prescribed(order);
cons_val_phi             =  prescribed;

