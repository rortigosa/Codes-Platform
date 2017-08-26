function  str   =  InternalWorkAssemblyFormulation(str)

switch str.data.formulation
    case 'u'
         str    =  ParallelInternalWorkUAssembly(str);
    case 'up'
         str    =  ParallelInternalWorkUPAssembly(str);
    case 'FHJ'
         str    =  ParallelInternalWorkFHJAssembly(str);
    case 'CGC' 
         str    =  ParallelInternalWorkCGCAssembly(str);
    case 'CGCCascade'
         str    =  ParallelInternalWorkCGCCascadeAssembly(str);
    case {'electro_mechanics','electro_mechanics_BEM_FEM'}
         str    =  ParallelInternalWorkElectroAssembly(str);
    case {'electro_mechanics_Helmholtz','electro_mechanics_BEM_FEM_Helmholtz'}
         str    =  ParallelInternalWorkElectroHelmholtzAssembly(str);
    case {'electro_mechanics_incompressible','electro_mechanics_incompressible_BEM_FEM'}
         str    =  ParallelInternalWorkElectroIncompressibleAssembly(str);
    case {'electro_mechanics_Helmholtz_incompressible','electro_mechanics_Helmholtz_incompressible_BEM_FEM'}
         str    =  ParallelInternalWorkElectroHelmholtzIncompressibleAssembly(str);
    case {'electro_mechanics_mixed_incompressible','electro_mechanics_mixed_incompressible_BEM_FEM'}        
         str    =  ParallelInternalWorkElectroMixedIncompressibleAssembly(str);
    case {'electro','electro_BEM_FEM'}        
         str    =  ParallelInternalWorkOnlyElectroAssembly(str);
end
