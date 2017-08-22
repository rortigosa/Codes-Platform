function [str]                                 =  postprocessing_acoustic_tensor(str,n_alpha,n_beta)

%--------------------------------------------------------------------------
% 1. Finite element technology and Gauss integration for postprocessing
%--------------------------------------------------------------------------
% if ~isfield(str,'postproc')
 str.postproc.output_type                      =  'averaging';
% end
%--------------------------------------------------------------------------
% 1.2 Reordering of shape functions and connectivities.
%--------------------------------------------------------------------------
switch str.data.dim
    case {13,23}
        str.Eulerian_x                         =  str.solid.BEAM_SHELL.continuum.Eulerian_x;
end
switch str.data.dim
    case {2,3}
         old_str.data.dim                      =  str.data.dim;
    case {13,23}
         old_str.data.dim                      =  str.data.dim;
         str.data.dim                          =  3;         
end

%--------------------------------------------------------------------------
% 1.2 Shape functions in the solution nodes of the elements.
%--------------------------------------------------------------------------
switch old_str.data.dim
    case {2,3}
         str.f_e                               =  str.solution_nodes.f_e;
         str.quadrature.Chi                    =  str.solution_nodes.quadrature.Chi;
    case {13,23}
         str.f_e                               =  str.postproc.f_e;
         str.quadrature                        =  str.solid.BEAM_SHELL.continuum.f_e.isoparametric_nodal_positions;
end
%-----------------------------------------------------------------
% 2. Initialisation of variables to obtain in the postprocessing.
%-----------------------------------------------------------------
switch str.postproc.output_type
    case 'averaging'
        str                          =  postpocessing_acoustic_tensor_initialisation(str);
    case 'discontinuous'
end
%-----------------------------------------------------------------
% 3. Postprocessing exclusively solutions in the postprocessing mesh.
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% 3. Postprocessing per element of all the postprocessing
%    variables in the postprocessing mesh.
%-----------------------------------------------------------------
str                                   =  postpocessing_computing_acoustic_tensor(str,n_alpha,n_beta);
%-----------------------------------------------------------------
%
%-----------------------------------------------------------------
[str]                                 =  plot_variables_treatment(str);
%-----------------------------------------------------------------
%  Express results in the postprocessing mesh.
%-----------------------------------------------------------------
%-----------------------------------------------------------------
% 6. Method to compute the variables in the vertices(output type): averaging...
%-----------------------------------------------------------------
[str]                                 =  average_acoustic_tensor(str);

str.data.dim                                   =  old_str.data.dim;