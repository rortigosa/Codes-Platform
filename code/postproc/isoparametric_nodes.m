%**************************************************************************
%**************************************************************************
%  1. This function needs the following information:
%     -str.data.shape             == 1 or 0 depending on if the element is 
%                                    quadrilateral or triangular.
%     -str.data.degree            == Degree of the polynomial used in the 
%                                    finite element interpolation.
%
%  2. This function will:
%  2.1 Store the Isoparametric coordinates of the Gauss points.
%  2.2 Store the weights used in the Gauss integration points.
%  In the case of 4 Gauss points, the information would be stored as
%  follows,
%  W  =  W_x_1*W_y_1      Points = chi_1    eta_1
%        W_x_2*W_y_1               chi_2    eta_2
%        W_x_1*W_y_2               chi_1    eta_2
%        W_x_2*W_y_2               chi_2    eta_2
%  
%**************************************************************************
%**************************************************************************


function [str]                       =  isoparametric_nodes(str)

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 1. Data needed
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
shape                                =  str.data.shape;

switch str.data.dim
    case 2
         %--------------------------------------------------------------------------
         %  2. Quadrilateral elements.
         %--------------------------------------------------------------------------
         if shape==1
             %--------------------------------------------------------------------------
             %  2.2.1 Points in the desired format
             %--------------------------------------------------------------------------
             chi_v                                =  [-1 1 1 -1];
             eta_v                                =  [-1 -1 1 1];
             str.quadrature.Chi                   =  [chi_v;eta_v]';
             %--------------------------------------------------------------------------
             %  3. Triangular elements.
             %--------------------------------------------------------------------------
         elseif shape==0
         end
   
    case 3
         if shape==1
             %--------------------------------------------------------------------------
             %  2.2.1 Points in the desired format
             %--------------------------------------------------------------------------
             chi_v                                =  [-1  1  1 -1  -1  1  1  -1 ];
             eta_v                                =  [-1 -1  1  1  -1 -1  1   1];
             iota_v                               =  [-1 -1 -1 -1   1  1  1   1];
             str.quadrature.Chi                   =  [chi_v;eta_v;iota_v]';
         end
end