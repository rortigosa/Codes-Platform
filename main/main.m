%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main program that will execute all the functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main
  
clc                              
warning('off','all')               
clear all
%load('prep1.mat')
%--------------------------------------------------------------------------
% Select jobfolder.
%--------------------------------------------------------------------------
%basedir_fem     =  'C:\MECHANICAL_SHELL_CODE_OPTIMISED\';
%basedir_fem     =  'C:\SoftwareDevelopment\NEW_ELECTROMECHANICS_CODE_v2\';
basedir_fem      =  'C:\SoftwareDevelopment\Codes-Platform\';

addpath(genpath(fullfile(basedir_fem,'code')));
addpath(genpath(fullfile(basedir_fem,'main')));
%--------------------------------------------------------------------------
% Select example
%--------------------------------------------------------------------------
input            =  ExampleSelection(basedir_fem);
str.jobfolder    =  input.jobfolder;
%--------------------------------------------------------------------------
% New preprocessor or load former one 
%--------------------------------------------------------------------------
newpreprocessor  =  1;
switch newpreprocessor
    case 1
         %-----------------------------------------------------------------
         % Add to path 
         %-----------------------------------------------------------------
         str     =  AddPathFunction(str,input);
         %-----------------------------------------------------------------
         % Initial data     
         %-----------------------------------------------------------------
         str     =  InitialData(str);
         %-----------------------------------------------------------------
         % Geometry preprocessor    
         %-----------------------------------------------------------------
         str     =  GeometryPreprocessor(str); 
         %-----------------------------------------------------------------
         % Gauss quadrature across each of the layers 
         %-----------------------------------------------------------------
         str     =  GetQuadratureRules(str);  
         %-----------------------------------------------------------------
         % Finite Element shape functions 
         %-----------------------------------------------------------------
         str     =  FEMShapeFunctions(str);
         %-----------------------------------------------------------------
         % Model geometrical information.   
         %-----------------------------------------------------------------
         str     =  MeshGenerationFormulation(str);
         %-----------------------------------------------------------------
         % Boundary preprocessing.     
         %-----------------------------------------------------------------
         str     =  BoundaryPreprocessingFinal(str);
         %-----------------------------------------------------------------
         % Material information.      
         %-----------------------------------------------------------------
         str     =  MaterialModelPreprocessor(str);
         %-----------------------------------------------------------------
         % Initialisation                  
         %-----------------------------------------------------------------
         str     =  InitialisationFormulation(str);
         %-----------------------------------------------------------------
         % Constraints                      
         %-----------------------------------------------------------------
         str     =  BoundaryConditionsManager(str);
    case 0
load('saved_preprocessor.mat');     
end  
%--------------------------------------------------------------------------
% Solver                                                                                                              
%-------------------------------------------------------------------------- 
Solver(str);
 
end

  


