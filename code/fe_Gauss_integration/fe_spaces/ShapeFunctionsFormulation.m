%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Types of shape functions: volume and surface
% Within, we can find: 
%
%  volume---->  bilinear, mass
%  surface--->  bilinear, mass, contact
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function  str                                =  ShapeFunctionsFormulation(str)

switch str.data.formulation
    case 'u'
          str                                =  ShapeFunctionsU(str);
    case 'up'
          str                                =  ShapeFunctionsUP(str);
    case 'FHJ'
          str                                =  ShapeFunctionsFHJ(str);
    case {'CGC','CGCCascade'}
          str                                =  ShapeFunctionsCGC(str);        
    case {'electro_mechanics','electro_mechanics_BEM_FEM','electro_mechanics_Helmholtz',...
          'electro_mechanics_Helmholtz_BEM_FEM','electro_BEM_FEM','electro'}
          str                                =  ShapeFunctionsElectro(str);
    case {'electro_mechanics_incompressible','electro_mechanics_incompressible_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
          str                                =  ShapeFunctionsElectroIncompressible(str);        
    case {'electro_mechanics_mixed_incompressible','electro_mechanics_mixed_incompressible_BEM_FEM'}
          str                                =  ShapeFunctionsElectroMixedIncompressible(str);                
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Shape functions for contact (2D shape functions).
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[N,DN_chi]                                   =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.contact,str.geometry.dim-1,str.quadrature.surface.contact);
str.fem.surface.contact.contact.N            =  N;
str.fem.surface.contact.contact.DN_chi       =  DN_chi;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Shape functions for electric flux in BEM/FEM problems
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
switch str.data.formulation
    case {'electro_BEM_FEM'}
         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.phi,str.geometry.dim-1,str.quadrature.surface.BEM_FEM);
         str.fem.surface.BEM_FEM.phi.N       =  N;
         str.fem.surface.BEM_FEM.phi.DN_chi  =  DN_chi;

         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.q,str.geometry.dim-1,str.quadrature.surface.BEM_FEM);
         str.fem.surface.BEM_FEM.q.N         =  N;
         str.fem.surface.BEM_FEM.q.DN_chi    =  DN_chi;

         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.q,str.geometry.dim-1,str.quadrature.surface.bilinear);
         str.fem.surface.bilinear.q.N        =  N;
         str.fem.surface.bilinear.q.DN_chi   =  DN_chi;
    case {'electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM','electro_mechanics_mixed_incompressible_BEM_FEM'}
         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.x,str.geometry.dim-1,str.quadrature.surface.BEM_FEM);
         str.fem.surface.BEM_FEM.x.N         =  N;
         str.fem.surface.BEM_FEM.x.DN_chi    =  DN_chi;
        
         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.phi,str.geometry.dim-1,str.quadrature.surface.BEM_FEM);
         str.fem.surface.BEM_FEM.phi.N       =  N;
         str.fem.surface.BEM_FEM.phi.DN_chi  =  DN_chi;

         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.q,str.geometry.dim-1,str.quadrature.surface.BEM_FEM);
         str.fem.surface.BEM_FEM.q.N         =  N;
         str.fem.surface.BEM_FEM.q.DN_chi    =  DN_chi;

         [N,DN_chi]                          =  ShapeFunctionComputation(str.fem.shape,str.fem.degree.q,str.geometry.dim-1,str.quadrature.surface.bilinear);
         str.fem.surface.bilinear.q.N        =  N;
         str.fem.surface.bilinear.q.DN_chi   =  DN_chi;
end
