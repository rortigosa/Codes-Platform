%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Compute second derivative of the internal energy with respect
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigenvalues                     =  eigenvalues_constitutive_tensor(str,F,H,J)

%F=eye(3);

H_vector                                 =  reshape(H',9,1);
E                                        =  zeros(3,3,3);
E(1,2,3)                                 =  1;
E(3,1,2)                                 =  1;
E(2,3,1)                                 =  1;
E(1,3,2)                                 =  -1;
E(3,2,1)                                 =  -1;
E(2,1,3)                                 =  -1;

%--------------------------------------------------------------------------    
% Derivatives for the case of a Neo-Hookean model
%--------------------------------------------------------------------------    
alpha                                    =  str.properties.material_parameters.alpha(1,:);
lambda                                   =  str.properties.material_parameters.lambda(1,:);
D2WDFDF                                  =  2*alpha*eye(9);
D2WDJDJ                                  =  2*alpha/J^2 + lambda;
DWDJ                                     =  -2*alpha/J + lambda*(J - 1);

%--------------------------------------------------------------------------    
% Geometric term
%--------------------------------------------------------------------------    
tensor                                   =  zeros(3,3,3,3);
for l=1:3
    for k=1:3
        for j=1:3
            for i=1:3
                for p=1:3
                    for q=1:3
                        tensor(i,j,k,l)  =  tensor(i,j,k,l) + E(i,k,p)*F(p,q)*E(q,j,l);
                    end
                end
            end
        end
    end
end
%tensor                                    =  permute(tensor,[1 3 2 4]);


%--------------------------------------------------------------------------
% Transformation of the nonlinear term into a matrix
%--------------------------------------------------------------------------
matrix(1,1)                              =  tensor(1,1,1,1);
matrix(1,2)                              =  tensor(1,1,1,2);
matrix(1,3)                              =  tensor(1,1,1,3);
matrix(1,4)                              =  tensor(1,1,2,1);
matrix(1,5)                              =  tensor(1,1,2,2);
matrix(1,6)                              =  tensor(1,1,2,3);
matrix(1,7)                              =  tensor(1,1,3,1);
matrix(1,8)                              =  tensor(1,1,3,2);
matrix(1,9)                              =  tensor(1,1,3,3);

matrix(2,1)                              =  tensor(1,2,1,1);
matrix(2,2)                              =  tensor(1,2,1,2);
matrix(2,3)                              =  tensor(1,2,1,3);
matrix(2,4)                              =  tensor(1,2,2,1);
matrix(2,5)                              =  tensor(1,2,2,2);
matrix(2,6)                              =  tensor(1,2,2,3);
matrix(2,7)                              =  tensor(1,2,3,1);
matrix(2,8)                              =  tensor(1,2,3,2);
matrix(2,9)                              =  tensor(1,2,3,3);

matrix(3,1)                              =  tensor(1,3,1,1);
matrix(3,2)                              =  tensor(1,3,1,2);
matrix(3,3)                              =  tensor(1,3,1,3);
matrix(3,4)                              =  tensor(1,3,2,1);
matrix(3,5)                              =  tensor(1,3,2,2);
matrix(3,6)                              =  tensor(1,3,2,3);
matrix(3,7)                              =  tensor(1,3,3,1);
matrix(3,8)                              =  tensor(1,3,3,2);
matrix(3,9)                              =  tensor(1,3,3,3);

matrix(4,1)                              =  tensor(2,1,1,1);
matrix(4,2)                              =  tensor(2,1,1,2);
matrix(4,3)                              =  tensor(2,1,1,3);
matrix(4,4)                              =  tensor(2,1,2,1);
matrix(4,5)                              =  tensor(2,1,2,2);
matrix(4,6)                              =  tensor(2,1,2,3);
matrix(4,7)                              =  tensor(2,1,3,1);
matrix(4,8)                              =  tensor(2,1,3,2);
matrix(4,9)                              =  tensor(2,1,3,3);

matrix(5,1)                              =  tensor(2,2,1,1);
matrix(5,2)                              =  tensor(2,2,1,2);
matrix(5,3)                              =  tensor(2,2,1,3);
matrix(5,4)                              =  tensor(2,2,2,1);
matrix(5,5)                              =  tensor(2,2,2,2);
matrix(5,6)                              =  tensor(2,2,2,3);
matrix(5,7)                              =  tensor(2,2,3,1);
matrix(5,8)                              =  tensor(2,2,3,2);
matrix(5,9)                              =  tensor(2,2,3,3);

matrix(6,1)                              =  tensor(2,3,1,1);
matrix(6,2)                              =  tensor(2,3,1,2);
matrix(6,3)                              =  tensor(2,3,1,3);
matrix(6,4)                              =  tensor(2,3,2,1);
matrix(6,5)                              =  tensor(2,3,2,2);
matrix(6,6)                              =  tensor(2,3,2,3);
matrix(6,7)                              =  tensor(2,3,3,1);
matrix(6,8)                              =  tensor(2,3,3,2);
matrix(6,9)                              =  tensor(2,3,3,3);

matrix(7,1)                              =  tensor(3,1,1,1);
matrix(7,2)                              =  tensor(3,1,1,2);
matrix(7,3)                              =  tensor(3,1,1,3);
matrix(7,4)                              =  tensor(3,1,2,1);
matrix(7,5)                              =  tensor(3,1,2,2);
matrix(7,6)                              =  tensor(3,1,2,3);
matrix(7,7)                              =  tensor(3,1,3,1);
matrix(7,8)                              =  tensor(3,1,3,2);
matrix(7,9)                              =  tensor(3,1,3,3);

matrix(8,1)                              =  tensor(3,2,1,1);
matrix(8,2)                              =  tensor(3,2,1,2);
matrix(8,3)                              =  tensor(3,2,1,3);
matrix(8,4)                              =  tensor(3,2,2,1);
matrix(8,5)                              =  tensor(3,2,2,2);
matrix(8,6)                              =  tensor(3,2,2,3);
matrix(8,7)                              =  tensor(3,2,3,1);
matrix(8,8)                              =  tensor(3,2,3,2);
matrix(8,9)                              =  tensor(3,2,3,3);

matrix(9,1)                              =  tensor(3,3,1,1);
matrix(9,2)                              =  tensor(3,3,1,2);
matrix(9,3)                              =  tensor(3,3,1,3);
matrix(9,4)                              =  tensor(3,3,2,1);
matrix(9,5)                              =  tensor(3,3,2,2);
matrix(9,6)                              =  tensor(3,3,2,3);
matrix(9,7)                              =  tensor(3,3,3,1);
matrix(9,8)                              =  tensor(3,3,3,2);
matrix(9,9)                              =  tensor(3,3,3,3);

%--------------------------------------------------------------------------
% Total derivative. 
%--------------------------------------------------------------------------
constitutive_tensor                      =  D2WDFDF + D2WDJDJ*(H_vector*H_vector') + DWDJ*matrix;

%--------------------------------------------------------------------------
% Eigenvalues of the matrix.
%--------------------------------------------------------------------------
eigenvalues                              =  sort(eig(constitutive_tensor));

