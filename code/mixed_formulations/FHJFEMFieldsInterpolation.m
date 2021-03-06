%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 
% In this function, FEM interpolation of the following fields is carried
% out: (F,H,J,D0,d), (SigmaF,SigmaH,SigmaJ,Sigmad) and pressure
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [F,Fvectorised,H,...
     Hvectorised,J,...
    SigmaF,SigmaFvectorised,...
    SigmaH,SigmaHvectorised,...
    SigmaJ]           =  FHJFEMFieldsInterpolation(ielem,dim,ngauss,fem,solution,mesh_F,mesh_H,mesh_J)

%--------------------------------------------------------------------------
% Interpolation of F field
%--------------------------------------------------------------------------
Felem_vectorised        =  solution.F(:,mesh_F.connectivity(:,ielem));
Fvectorised             =  VectorFEMInterpolation(ngauss,fem.volume.bilinear.F.N,Felem_vectorised);
F                       =  permute((reshape(Fvectorised,dim,dim,[])),[2 1 3]);
%--------------------------------------------------------------------------
% Interpolation of H field
%--------------------------------------------------------------------------
Helem_vectorised        =  solution.H(:,mesh_H.connectivity(:,ielem));
Hvectorised             =  VectorFEMInterpolation(ngauss,fem.volume.bilinear.H.N,Helem_vectorised);
H                       =  permute((reshape(Hvectorised,dim,dim,[])),[2 1 3]);
%--------------------------------------------------------------------------
% Interpolation of J field
%--------------------------------------------------------------------------
J                       =  ScalarFEMInterpolation(fem.volume.bilinear.J.N,...
                                                 solution.J(mesh_J.connectivity(:,ielem)));                                         
%--------------------------------------------------------------------------
% Interpolation of SigmaF field
%--------------------------------------------------------------------------                                             
SigmaFelem_vectorised   =  solution.SigmaF(:,mesh_F.connectivity(:,ielem));
SigmaFvectorised        =  VectorFEMInterpolation(ngauss,fem.volume.bilinear.F.N,SigmaFelem_vectorised);
SigmaF                  =  permute((reshape(SigmaFvectorised,dim,dim,[])),[2 1 3]);
%--------------------------------------------------------------------------
% Interpolation of SigmaH field
%--------------------------------------------------------------------------                                             
SigmaHelem_vectorised   =  solution.SigmaH(:,mesh_H.connectivity(:,ielem));
SigmaHvectorised        =  VectorFEMInterpolation(ngauss,fem.volume.bilinear.H.N,SigmaHelem_vectorised);
SigmaH                  =  permute((reshape(SigmaHvectorised,dim,dim,[])),[2 1 3]);
%--------------------------------------------------------------------------
% Interpolation of SigmaJ field
%--------------------------------------------------------------------------                                             
SigmaJ                  =  ScalarFEMInterpolation(fem.volume.bilinear.J.N,...
                                                 solution.SigmaJ(mesh_J.connectivity(:,ielem)));                                         
