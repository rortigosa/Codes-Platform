%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% We specify nodal loads
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function  P_nodal  =  NodalLoads(str)

zmin               =  min(str.mesh.volume.x.nodes(3,:));
X0                 =  [0;0;zmin];
for inode=1:str.mesh.volume.x.n_nodes
    if norm(str.mesh.volume.x.nodes(:,inode)-X0)<1e-6
       node        =  inode;
       break;
    end 
end

fixdofs            =  3*node;
consval            =  0*5e3*ones(size(fixdofs,1),1);

n_dofs             =  size(str.solution.x.Eulerian_x(:),1);
P_nodal            =  zeros(n_dofs,1);
P_nodal(fixdofs)   =  consval;


