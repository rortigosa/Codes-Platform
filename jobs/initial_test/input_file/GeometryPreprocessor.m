

function str                =  GeometryPreprocessor(str)

str.geometry.type           =  'Structured_Extruded_Hexa_Prism';

switch str.geometry.type
    case 'Structured_Triad_Rectangle'    
          str.geometry.StructuredTriadRectangle.Lx   =  1;
          str.geometry.StructuredTriadRectangle.Ly   =  1;
          str.geometry.StructuredTriadRectangle.Nx   =  10;
          str.geometry.StructuredTriadRectangle.Ny   =  10;
          str.geometry.dim                           =  2;
    case 'Structured_Quad_Rectangle'    
          str.geometry.StructuredQuadRectangle.Lx   =  1;
          str.geometry.StructuredQuadRectangle.Ly   =  1;
          str.geometry.StructuredQuadRectangle.Nx   =  10;
          str.geometry.StructuredQuadRectangle.Ny   =  10;
          str.geometry.dim                          =  2;
    case 'Structured_Extruded_Hexa_Prism'
          str.geometry.StructuredExtrudedHexaPrism.Lx          =  1;
          str.geometry.StructuredExtrudedHexaPrism.Ly          =  1;
          str.geometry.StructuredExtrudedHexaPrism.Nx          =  4;
          str.geometry.StructuredExtrudedHexaPrism.Ny          =  4;
          str.geometry.StructuredExtrudedHexaPrism.Nz          =  4;
          thickness                                            =  1;
          str.geometry.StructuredExtrudedHexaPrism.thickness   =  repmat(thickness/str.geometry.StructuredExtrudedHexaPrism.Nz,str.geometry.StructuredExtrudedHexaPrism.Nz,1);
          str.geometry.dim                                     =  3;
    case 'Structured_Tet_Prism'
          str.geometry.StructuredTetPrism.Lx   =  1;
          str.geometry.StructuredTetPrism.Ly   =  1;
          str.geometry.StructuredTetPrism.Lz   =  1;
          str.geometry.StructuredTetPrism.Nx   =  10;
          str.geometry.StructuredTetPrism.Ny   =  10;
          str.geometry.StructuredTetPrism.Nz   =  10;
          str.geometry.dim                     =  3;
    case 'Structured_Hexa_Prism'
          str.geometry.StructuredTetPrism.Lx   =  1;
          str.geometry.StructuredTetPrism.Ly   =  1;
          str.geometry.StructuredTetPrism.Lz   =  1;
          str.geometry.StructuredTetPrism.Nx   =  10;
          str.geometry.StructuredTetPrism.Ny   =  10;
          str.geometry.StructuredTetPrism.Nz   =  10;
          str.geometry.dim                     =  3;
    case 'Structured_Quad_Circle'
          str.geometry.StructuredQuadCircle.r         =  1;
          str.geometry.StructuredQuadCircle.Nr        =  1;
          str.geometry.StructuredQuadCircle.Ntheta    =  1;
          str.geometry.dim                            =  2;
    case 'Structured_Extruded_Hexa_Cilinder'
          str.geometry.StructuredExtrudedHexaCilinder.r          =  1;
          str.geometry.StructuredExtrudedHexaCilinder.Nr         =  10;
          str.geometry.StructuredExtrudedHexaCilinder.Ntheta     =  10;
          str.geometry.StructuredExtrudedHexaCilinder.Nz         =  10;
          thickness                                              =  1;
          str.geometry.StructuredExtrudedHexaCilinder.thickness  =  repmat(thickness/str.geometry.StructuredExtrudedHexaPrism.Nz,str.geometry.StructuredExtrudedHexaPrism.Nz,1);
          str.geometry.dim                                       =  3;
end        
