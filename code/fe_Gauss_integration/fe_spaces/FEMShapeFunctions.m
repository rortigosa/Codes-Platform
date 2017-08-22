%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Create shape functions
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str    =  FEMShapeFunctions(str)

str             =  ShapeFunctionsFormulation(str);
str             =  ShapeFunctionNodes(str);
str             =  ShapeFunctionPostprocessingNodes(str);
