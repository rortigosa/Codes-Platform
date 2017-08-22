function [J,J0]       =  Intensity_vector(str,D0,D0n1,Jn,Fn)

J0                    =  zeros(str.data.dim,size(str.quadrature.Chi,1));
J                     =  J0;
Dt                    =  str.data.delta_t;
for inode=1:size(str.quadrature.Chi,1)
    %---------------------------------------------------------------------
    % Obtain Lagrangian electric displacement.
    %---------------------------------------------------------------------
    D_                =  1/Jn(inode,1)*Fn(:,:,inode)*D0(:,inode);
    Dn1_              =  1/Jn(inode,1)*Fn(:,:,inode)*D0n1(:,inode);
    DDt               =  (D_ - Dn1_)/Dt;
    %-------------------------------------------------------------------
    % Density vector.
    %-------------------------------------------------------------------
    J(:,inode)          =  -DDt;
    J0(:,inode)         =  Jn(inode,1)*Fn(:,:,inode)\J(:,inode);
end

