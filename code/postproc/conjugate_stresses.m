%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cauchy stress for postprocessing. We have to be very carefull because the
% resulting Cauchy stress needs to satisfy material frame indefference.
% Therefore, it needs to be symmetric.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SigmaF, SigmaH, SigmaJ]                                  = conjugate_stresses(str,F_,H_,J_)

SigmaF                                =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
SigmaH                                =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
SigmaJ                                =  zeros(size(str.postproc.f_e.N,2),1);


for igauss=1:size(str.postproc.f_e.N,2)
    F                                =  F_(:,:,igauss);
    H                                =  H_(:,:,igauss);
    J                                =  J_(igauss);
switch str.material_model_info.material_model
    case {'Saint-Venant'}
         lambda                      =  str.properties.lambda(str.properties.material_identifier);
         mu                          =  str.properties.mu(str.properties.material_identifier);
         SigmaF(:,:,igauss)          =  mu*(trace(F'*F) - 1)*F + lambda/2*(trace(F'*F) - 3);
         SigmaH(:,:,igauss)          =  -mu*H;
    case {'Neo-Hookean'}
         lambda                      =  str.properties.lambda(str.properties.material_identifier);
         mu                          =  str.properties.mu(str.properties.material_identifier);
         SigmaF(:,:,igauss)          =  mu*F;
         SigmaJ(igauss,1)            =  -mu/J + lambda*(J - 1);
    case 'Mooney-Rivlin'
         alpha                       =  str.properties.alpha(str.properties.material_identifier);
         beta                        =  str.properties.beta(str.properties.material_identifier);
         lambda                      =  str.properties.lambda(str.properties.material_identifier);

         SigmaF(:,:,igauss)          =  alpha*F;
         SigmaH(:,:,igauss)          =  beta*H;
         SigmaJ(igauss,1)            =  -(alpha + 2*beta)/J + lambda*(J - 1);
    case 'Modified-Mooney-Rivlin'        
         alpha                       =  str.properties.alpha(str.properties.material_identifier);
         beta                        =  str.properties.beta(str.properties.material_identifier);
         lambda                      =  str.properties.lambda(str.properties.material_identifier);
        
         SigmaF(:,:,igauss)          =  alpha*trace(F'*F)*F;
         SigmaH(:,:,igauss)          =  beta*trace(H'*H)*H;
         SigmaJ(igauss,1)            =  -(alpha + 2*beta)/J + lambda*(J - 1);

end
end
