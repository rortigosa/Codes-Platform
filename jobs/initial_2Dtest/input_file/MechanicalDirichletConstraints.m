%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  In this function we specify boundary conditions for the mechanical
%  physics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    [fixdof_u,freedof_u,...
    cons_val_u]         =  MechanicalDirichletConstraints(str)

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
nodes1                  =  zeros(str.mesh.volume.x.n_nodes,1);
presc1                  =  zeros(str.mesh.volume.x.n_nodes,1);
xmin                    =  min(str.mesh.volume.x.nodes(1,:));
for inode=1:str.mesh.volume.x.n_nodes
    x                   =  str.mesh.volume.x.nodes(:,inode);
    if (x(1) - xmin)<1e-6
       nodes1(inode)    =  inode;
       presc1(inode)    =  0;
    end 
end
nodes1                  =  nodes1(nodes1>0);
presc1                  =  presc1(nodes1);
%--------------------------------------------------------------------------
% Degrees of freedom with associated constraints
%--------------------------------------------------------------------------
dofs1x                  =  1 + 3*(nodes1 - 1);  %  Select direction X for set 1
dofs1y                  =  2 + 3*(nodes1 - 1);  %  Select direction Y for set 1

freedof_u               =  (1:2*str.mesh.volume.x.n_nodes)';
[fixdof_u,order]        =  sort([dofs1x;dofs1y]);
freedof_u(fixdof_u)     =  [];
%--------------------------------------------------------------------------
% Degrees of freedom with associated constraints 
%--------------------------------------------------------------------------
presc1x                 =  presc1;
presc1y                 =  presc1;

prescribed              =  [presc1x;presc1y];
prescribed              =  prescribed(order);
cons_val_u              =  prescribed;

