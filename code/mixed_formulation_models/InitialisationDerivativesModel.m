function str  =  InitialisationDerivativesModel(str)

switch str.data.formulation
    case {'u','up','FHJ'}
       str    =  InitialisationDerivativesMechanics(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);
    case {'CGC','CGCCascade'}
       str    =  InitialisationDerivativesMechanicsC(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);
    case {'electro_mechanics','electro_mechanics_incompressible','electro_mechanics_BEM_FEM','electro_mechanics_incompressible_BEM_FEM'}
       str    =  InitialisationDerivativesElectro(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);      
    case {'electro_mechanics_mixed_incompressible','electro_mechanics_mixed_incompressible_BEM_FEM'}
       str    =  InitialisationDerivativesElectroMixed(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);      
    case {'electro_mechanics_Helmholtz','electro_mechanics_Helmholtz_BEM_FEM',...
          'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
       str    =  InitialisationDerivativesElectroHelmholtz(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);      
    case {'electro_BEM_FEM'}
       str    =  InitialisationDerivativesOnlyElectroBEMFEM(str.geometry.dim,size(str.quadrature.volume.bilinear.Chi,1),str);      
end