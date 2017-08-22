function [nodes_3D,...
          connectivity_3D]     =  StructuredExtrudedHexaPrismMeshGenerator(degree,Lx,Ly,...
                                                                 Nx,Ny,Nz,thickness,Nshape)

[nodes,connectivity]           =  StructuredQuadRectangleMeshGenerator(degree,Lx,Ly,Nx,Ny);
[nodes_3D,connectivity_3D]     =  ThreeDSolidShellReconstruction(nodes',connectivity',degree,Nz,thickness,Nshape);


