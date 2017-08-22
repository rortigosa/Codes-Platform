function [str]                                     =  average_acoustic_tensor(str)


nodes_counter                                      =  str.nodes_counter;
switch str.postproc.output_type
    case 'averaging' 
         str.postproc.q               =  str.postproc.q'./nodes_counter;
end