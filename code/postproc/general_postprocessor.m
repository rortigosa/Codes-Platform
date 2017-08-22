function [str]                                   =  general_postprocessor(str)


jobfolder                                        =  (([str.jobfolder '\results']));
cd(jobfolder);
filename                                         =  ['time' num2str(str.time) '.mat'];
if str.time_iteration>1400
save(filename);
end
str.postproc.displacement(str.time_iteration)     =  str.Eulerian_x(1,35) - str.Lagrangian_X(1,35);

%----------------------------
% Tecplot postprocessing.
%----------------------------
%TECPLOT_plotting(str)
%----------------------------
% Intensity and electric power.
%----------------------------
if 0
if isfield(str.postproc, 'n_Intensity_regions')
    if str.time_iteration>=2
        [I,W]                                    =  Electric_Intensity(str);
        [I_analytical,W_analytical]              =  analytical_Electric_Intensity(str);
        str.I(str.time_iteration)                =  I;
        str.I_analytical(str.time_iteration)     =  I_analytical;
        str.W(str.time_iteration)                =  W;
        str.W_analytical(str.time_iteration)     =  W_analytical;
        figure(300)
        plot(str.I(2:str.time_iteration),'-*')
        hold on
        plot(str.I_analytical(2:str.time_iteration),'r-*')
        figure(400)
        plot(str.W(2:str.time_iteration),'-*')
        hold on
        plot(str.W_analytical(2:str.time_iteration),'r-*')
    end
end
end
