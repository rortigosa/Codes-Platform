function [eul_D,E,E_norm] =  elec_disp_paraview(str,D0,E0)
    eul_D                 =  zeros(str.data.dim,size(str.quadrature.Chi,1));
    E                     =  zeros(str.data.dim,size(str.quadrature.Chi,1));
    E_norm                =  zeros(size(str.quadrature.Chi,1),1);
    for igauss=1:size(str.postproc.f_e.N,2)
        F                 =  str.grad.F(:,:,igauss);
        J                 =  det(F);
        eul_D(:,igauss)   =  1/J*F*D0(:,igauss);
        E(:,igauss)       =  inv(F)'*E0(:,igauss);  
        E_norm(igauss,1)  =  norm(E(:,igauss));
    end
end
