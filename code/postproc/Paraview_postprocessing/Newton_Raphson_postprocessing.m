%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Newton-Raphson postprocessing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str                     =  Newton_Raphson_postprocessing(str,incr_load,newton_raphson_iteration,acumulated_factor,load_factor)

figure(1)
plot(log10(str.Residual.Residual_stored(incr_load,1:newton_raphson_iteration)),'b-o')
switch str.data.analysis
    case 'static'
        str.temp.incr_load       =  incr_load;
        str                      =  saving_to_memory(str,acumulated_factor,incr_load,load_factor);
        jobfolder                =  (([str.jobfolder '\results']));
        cd(jobfolder);
        filename                 =  ['Load_increment_' num2str(incr_load)];
        save(filename,'-v7.3');
        str.time                 =  acumulated_factor;
        
        str.new_connectivity = str.connectivity;
        str.new_connectivity = str.new_connectivity(:,[1,10,3,6,7,8,2,4,9,5]);
%           str.new_connectivity = str.new_connectivity(:,[1 3 6 2 5 4]);
        str.mod_displacements = sqrt(str.Eulerian_x(:,1).^2+str.Eulerian_x(:,2).^2+str.Eulerian_x(:,3).^2);
        
        VTK_plotting('salida1.vtk','unstructured_grid',str.nodes,str.new_connectivity,'scalars','Desplazamiento',str.mod_displacements)
        VTK_plotting('salida2.vtk','unstructured_grid',str.Eulerian_x',str.new_connectivity,'scalars','Desplazamiento',str.mod_displacements)
        VTK_plotting('salidaP.vtk','unstructured_grid',str.nodes_p,str.connectivity_p,'scalars','Presion','Presionn',str.p,str.p)
%        TECPLOT_plotting(str);
end
