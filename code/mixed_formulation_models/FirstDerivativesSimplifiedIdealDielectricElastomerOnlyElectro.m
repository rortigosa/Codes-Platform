function mat_info  =  FirstDerivativesSimplifiedIdealDielectricElastomerOnlyElectro(ielem,...
                                                                                     D0,mat_info)
%--------------------------------------------------------------------------                                                                                        
% Material number identifier
%--------------------------------------------------------------------------                                                                                        
mat_id                    =  mat_info.material_identifier(ielem)                                                                                        ;
%--------------------------------------------------------------------------                                                                                        
% Material parameters
%--------------------------------------------------------------------------                                                                                        
e1                        =  mat_info.material_parameters.e2(mat_id);
%--------------------------------------------------------------------------                                                                                        
% First derivatives of the model
%--------------------------------------------------------------------------                                                                                        
mat_info.derivatives.DU.DUDD0    =  1/(e1)*D0;
