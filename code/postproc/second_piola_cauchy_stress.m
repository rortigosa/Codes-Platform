%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the Cauchy stress, we get Second Piola.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S]       =  second_piola_cauchy_stress(str,sigma)
S                  =  zeros(str.data.dim,str.data.dim,size(str.quadrature.Chi,1));
for iloop=1:size(str.quadrature.Chi,1)
    F              =  str.grad.F(:,:,iloop);
    invF           =  inv(F);
    J              =  det(F);
    S(:,:,iloop)   =  J*invF*sigma(:,:,iloop)*invF';
end
end
