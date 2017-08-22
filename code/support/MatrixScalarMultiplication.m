%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Compute the product between 
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function newmatrix       =  MatrixScalarMultiplication(dim,n_gauss,matrix,scalar)

newmatrix                =  zeros(dim,dim,n_gauss);
for igauss=1:n_gauss
    newmatrix(:,igauss)  =  matrix(:,:,igauss)*scalar(igauss);
end


