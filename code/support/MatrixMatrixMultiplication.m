%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Compute the product between 
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function newmatrix         =  MatrixMatrixMultiplication(dim,n_gauss,matrixA,matrixB)

newmatrix                  =  zeros(dim,dim,n_gauss);
for igauss=1:n_gauss
    newmatrix(:,:,igauss)  =  matrixA(:,:,igauss)*matrixB(:,:,igauss);
end


