%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Types of quadrature rules: volume and surface
% Within, we can find: 
%
%  volume---->  bilinear, mass
%  surface--->  bilinear, mass, contact
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                                  =  GetQuadratureRules(str) 

newstr                                        =  GaussQuadrature(str,str.quadrature.degree,str.geometry.dim-1);
str.quadrature.surface.bilinear.Chi           =  newstr.quadrature.Chi;
str.quadrature.surface.bilinear.W_v           =  newstr.quadrature.W_v; 

newstr.quadrature.Chi                         =  [];
newstr.quadrature.W_v                         =  [];
newstr                                        =  GaussQuadrature(newstr,str.quadrature.degree,str.geometry.dim);
str.quadrature.volume.bilinear.Chi            =  newstr.quadrature.Chi;
str.quadrature.volume.bilinear.W_v            =  newstr.quadrature.W_v;

newstr.quadrature.Chi                         =  [];
newstr.quadrature.W_v                         =  [];
newstr                                        =  GaussQuadrature(newstr,str.quadrature.degree*2,str.geometry.dim);
str.quadrature.volume.mass.Chi                =  newstr.quadrature.Chi;
str.quadrature.volume.mass.W_v                =  newstr.quadrature.W_v;

newstr.quadrature.Chi                         =  [];
newstr.quadrature.W_v                         =  [];
newstr                                        =  GaussQuadrature(newstr,str.quadrature.degree*2,str.geometry.dim-1);
str.quadrature.surface.mass.Chi               =  newstr.quadrature.Chi;
str.quadrature.surface.mass.W_v               =  newstr.quadrature.W_v;

newstr.quadrature.Chi                         =  [];
newstr.quadrature.W_v                         =  [];
newstr                                        =  GaussQuadrature(newstr,str.fem.degree.contact + str.fem.degree.x,str.geometry.dim-1);
str.quadrature.surface.contact.Chi            =  newstr.quadrature.Chi;
str.quadrature.surface.contact.W_v            =  newstr.quadrature.W_v; 

switch str.data.formulation
    case {'electro_BEM_FEM','electro_incompressible_BEM_FEM','electro_mixed_incompressible_BEM_FEM'}
         newstr.quadrature.Chi                =  [];
         newstr.quadrature.W_v                =  [];
         newstr                               =  GaussQuadrature(newstr,str.quadrature.degree,str.geometry.dim-1);
         str.quadrature.surface.BEM_FEM.Chi   =  newstr.quadrature.Chi;
         str.quadrature.surface.BEM_FEM.W_v   =  newstr.quadrature.W_v;
end
