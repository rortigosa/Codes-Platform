%**************************************************************************
%**************************************************************************
%  1.This function plots the electric potential distribution 
%**************************************************************************
%**************************************************************************

function  str=TECPLOT_plotting(str)   
if 1
    str.formulation                   =  str.old_formulation;
    str.old                           =  str;    
    switch str.vacuum.flag
        case 0
             % [str]                  =  flustr2plt_v2([],str,'str_tecplot.plt','fem_electromech',1,0);
             [str]                    =  provisional_flustr2plt_v2([],str,'TECPLOT_final.plt','fem_electromech',1,0);
             %save('prov')
        case 1 
             [str]                    =  provisional_flustr2plt_v2([],str,'str_tecplot_vacuum.plt','fem_electromech',1,0);
    end          
    switch str.data.post_storage 
           case 1
                %str                   =  str.postproc;    
        otherwise  
%    str          =  str.old;  % remove
    end
%     switch str.data.post_storage
%            case 1
%                 filename              =  ['postprocessing_results_time_' num2str(str.time) '.mat'];
%                 save(filename);
%     end
    
else
orderN                                =  [1 2 4 3];   
% %------------------------------------------
% % Plot electric potential
% %------------------------------------------
figure(2000)
for ielem=1:str.n_elem
        order(1:length(orderN))       =   str.connectivity(ielem,orderN);
        phielem                       =   str.phi(order);                
        xplot                         =   str.Eulerian_x(:,order);        
%        xplot                        =   str.Lagrangian_X(:,order);        
        hv=fill(xplot(1,:),xplot(2,:),phielem,'EdgeColor','none')  ;
%        hv=fill(xplot(1,:),xplot(2,:),phielem)  ;
    hold on
end  
colorbar('vert')
set(hv,'edgecolor','none'); 
set(hv,'MarkerEdgeColor','none');
set(hv); 
drawnow;
end








end
