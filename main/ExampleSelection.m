%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function selects the jobfolder of the example to analyse.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function input    =  ExampleSelection(basedir_fem)

%--------------------------------------------------------------------------
% Select the jogfolder 
%--------------------------------------------------------------------------
folder            =  'jobs\initial_test';
%--------------------------------------------------------------------------
% Input structure
%--------------------------------------------------------------------------
input.jobfolder   =  fullfile(basedir_fem,folder);


