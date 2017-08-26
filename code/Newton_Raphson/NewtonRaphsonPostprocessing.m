%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Newton-Raphson postprocessing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewtonRaphsonPostprocessing(str)

switch str.NR.convergence_plotting
    case 1
         figure(1)
         plot(log10(str.assembly.Residual_stored{str.NR.incr_load}(1:str.NR.iteration)),'b-o')
end
switch str.data.analysis
    case 'static'
        jobfolder                =  (([str.jobfolder '\results']));
        cd(jobfolder);
        filename                 =  ['Load_increment_' num2str(str.NR.incr_load)];
        save(filename,'-v7.3');
end
