%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  1.This function plots in TECPLOT 360 the deformed shape and the least of
%  the minors of the acoustic tensor
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function  str=TECPLOT_acoustic_tensor(str)   
  

[n_alpha,n_beta,Nodes_alpha,... 
    connectivity_alpha]               =  spherical_parametrisation_mesh;

%load(.....)       
          
        
  
%load('DG_mesh.mat')    
%str.postproc=alternative_str.postproc;
%str.nodes_counter = alternative_str.postproc.nodes_counter;          
  
[str]                                 =  acoustic_info([],str,n_alpha,n_beta,'acoustic_tensor_simulation.plt','fem_electromech',1,0);
   