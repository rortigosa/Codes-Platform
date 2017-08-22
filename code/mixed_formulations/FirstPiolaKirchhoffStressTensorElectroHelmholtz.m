%--------------------------------------------------------------------------
%--------------------------------------------------------------------------%
%   This function computes the first Piola-Kirchhoff stress tensor 
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function   Piola         =  FirstPiolaKirchhoffStressTensorElectroHelmholtz(ngauss,dim,F,H,DPsiDF,DPsiDH,DPsiDJ)

Piola                    =  zeros(dim,dim,ngauss);
for igauss=1:ngauss
    Piola_F              =  DPsiDF(:,:,igauss);
    Piola_H              =  JavierDoubleCrossProduct(DPsiDH(:,:,igauss),F(:,:,igauss));
    Piola_J              =  DPsiDJ(igauss)*H(:,:,igauss);
    Piola(:,:,igauss)    =  Piola_F + Piola_H + Piola_J;
end