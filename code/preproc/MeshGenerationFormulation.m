%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Read geometry and connectivities 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str                                =  MeshGenerationFormulation(str)

switch str.data.formulation
    case 'u'
         str                                =  MeshU(str);
    case 'up'
         str                                =  MeshUP(str);
    case 'FHJ'
         str                                =  MeshFHJ(str);
    case {'CGC','CGCCascade'}
         str                                =  MeshCGC(str);
    case {'electro_mechanics','electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz',...
          'electro_mechanics_Helmholtz_BEM_FEM','electro_BEM_FEM','electro'}
         str                                =  MeshElectro(str);
    case {'electro_mechanics_incompressible','electro_mechanics_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
         str                                =  MeshElectroIncompressible(str);        
    case {'electro_mechanics_mixed_incompressible','electro_mechanics_mixed_incompressible_BEM_FEM'}
         str                                =  MeshElectroMixedIncompressible(str);        
end
%--------------------------------------------------------------------------
% Contact
%--------------------------------------------------------------------------
if str.contact.lagrange_multiplier
    switch str.contact.method
        case 'mortar'
% [nodes,connectivity]                                  =  mesh_information_example(str.fem.mid_plane.contact.degree,'continuous',str);
% str.mesh.mid_plane.contact_multiplier.nodes           =  nodes';
% str.mesh.mid_plane.contact_multiplier.connectivity    =  connectivity';
% str.mesh.mid_plane.contact_multiplier.n_nodes         =  size(nodes,1);
% str.mesh.mid_plane.contact_multiplier.n_node_elem     =  size(connectivity,2);
        case 'node'
% str.mesh.mid_plane.contact_multiplier.nodes           =  str.mesh.mid_plane.x0.nodes;
% str.mesh.mid_plane.contact_multiplier.connectivity    =  str.mesh.mid_plane.x0.connectivity;
% str.mesh.mid_plane.contact_multiplier.n_nodes         =  str.mesh.mid_plane.x0.n_nodes;
% str.mesh.mid_plane.contact_multiplier.n_node_elem     =  str.mesh.mid_plane.x0.n_node_elem;            
    end
end
%--------------------------------------------------------------------------
% Mesh for BEM/FEM problems for the electric flux (continuous meshing for the moment)
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_BEM_FEM','electro_mechanics_Helmholtz_incompressible_BEM_FEM',...
          'electro_BEM_FEM'}
         [nodes,connectivity]              =  MeshGenerator(str,str.fem.degree.q,'continuous',str.fem.surface.nodes.q.DN_chi);         
         str.mesh.volume.q.nodes           =  nodes;
         str.mesh.volume.q.connectivity    =  connectivity;
         str.mesh.volume.q.n_nodes         =  size(nodes,2);
         str.mesh.volume.q.n_node_elem     =  size(connectivity,1);
end
