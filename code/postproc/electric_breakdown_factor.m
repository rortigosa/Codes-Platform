%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cauchy stress for postprocessing. We have to be very carefull because the
% resulting Cauchy stress needs to satisfy material frame indefference.
% Therefore, it needs to be symmetric.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function factor                                  = electric_breakdown_factor(str,E_norm,F_)


factor                          =  zeros(size(str.postproc.f_e.N,2),1);
if 0
    for igauss=1:size(str.postproc.f_e.N,2)
    switch str.data.formulation_type
        case {'displacement_potential_formulation','u_p_phi_compressible'}
            F                    =  str.grad.F(:,:,igauss);
        otherwise
        F                       =  F_(:,:,igauss);
        C                       =  F'*F;
    end
    lambda                      =  sort(sqrt((eig(C))));
    dimension                   =  max(size(lambda));
    lambda                      =  lambda(2:dimension);
    lambda                      =  sqrt(prod(lambda));
    if lambda<5.7
       E                        =  30.6*lambda^(1.13)*1e6;
    else
       E                        =  218*1e6;
    end
    factor(igauss)              =  E_norm(igauss)/E;
    end
end 
       
 
%asdf=98
    
