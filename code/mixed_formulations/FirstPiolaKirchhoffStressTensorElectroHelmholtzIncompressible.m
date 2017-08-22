%--------------------------------------------------------------------------
%--------------------------------------------------------------------------%
%   This function computes the first Piola-Kirchhoff stress tensor 
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function   Piola         =  FirstPiolaKirchhoffStressTensorElectroHelmholtzIncompressible(ngauss,dim,F,H,DPsiDF,DPsiDH,DPsiDJ,pressure)

Piola                    =  zeros(dim,dim,ngauss);
for igauss=1:ngauss
    Piola_F              =  DPsiDF(:,:,igauss);
    Piola_H              =  JavierDoubleCrossProduct(DPsiDH(:,:,igauss),F(:,:,igauss));
    Piola_J              =  DPsiDJ(igauss)*H(:,:,igauss);
    Piola_p              =  pressure(igauss)*H(:,:,igauss);
    Piola(:,:,igauss)    =  Piola_F + Piola_H + Piola_J + Piola_p;
end