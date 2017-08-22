%--------------------------------------------------------------------------
%--------------------------------------------------------------------------%
%   This function computes the first Piola-Kirchhoff stress tensor 
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function   Piola        =  FirstPiolaKirchhoffStressTensorU(ngauss,dim,F,H,DUDF,DUDH,DUDJ)

Piola                   =  zeros(dim,dim,ngauss);
for igauss=1:ngauss
    Piola_F             =  DUDF(:,:,igauss);
    Piola_H             =  JavierDoubleCrossProduct(DUDH(:,:,igauss),F(:,:,igauss));
    Piola_J             =  DUDJ(igauss)*H(:,:,igauss);
    Piola(:,:,igauss)   =  Piola_F + Piola_H + Piola_J;
end
