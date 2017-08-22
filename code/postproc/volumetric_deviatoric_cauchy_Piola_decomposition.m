function [sigma_dev,sigma_vol,...
    S_dev,S_vol,p]                  =  volumetric_deviatoric_cauchy_Piola_decomposition(str,sigma)

switch str.data.dim
    case 2
        dim                         =  2;
    case {3,13,23}
        dim                         = 3;
end
switch str.data.formulation_type
    %------------------------------------------------------------------
    %  Displacement potential formulation
    %------------------------------------------------------------------
    case {'displacement_potential_formulation','u_p_phi_compressible','twoD_u_V_mixed_formulation','twoD_u_V_v2_mixed_formulation','twoD_u_V_v4_mixed_formulation','twoD_u_V_v4_linearised_mixed_formulation'}
        sigma_dev               =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
        sigma_vol               =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
        S_dev                   =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
        S_vol                   =  zeros(str.data.dim,str.data.dim,size(str.postproc.f_e.N,2));
        p                       =  zeros(size(str.postproc.f_e.N,2),1);
        I                       =  eye(str.data.dim);
        for igauss=1:size(str.postproc.f_e.N,2)
        %-------------------------------------------------
        % sigma decomposition
        %-------------------------------------------------
        F                       =  str.grad.F(:,:,igauss);
        J                       =  str.grad.J(igauss);
        pressure                =  trace(sigma(:,:,igauss))/dim;
        sigma_vol_              =  pressure*I;
        sigma_dev_              =  sigma(:,:,igauss) - sigma_vol_;
        %-------------------------------------------------
        % S decomposition
        %-------------------------------------------------
        invF                    =  inv(F);
        S                       =  J*(F\(sigma(:,:,igauss)*inv(F)'));
        S_vol_                  =  J*pressure*invF*invF';
        S_dev_                  =  S - S_vol_;
        %-------------------------------------------------
        % Storage.
        %-------------------------------------------------
        sigma_vol(:,:,igauss)   =  sigma_vol_;
        sigma_dev(:,:,igauss)   =  sigma_dev_;
        p(igauss,1)             =  pressure;
        S_vol(:,:,igauss)       =  S_vol_;
        S_dev(:,:,igauss)       =  S_dev_;
        end
    case 'mixed_formulation'
        sigma_dev               =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
        sigma_vol               =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
        S_dev                   =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
        S_vol                   =  zeros(str.data.continuum_dim,str.data.continuum_dim,size(str.postproc.f_e.N,2));
        p                       =  zeros(size(str.postproc.f_e.N,2),1);
        I                       =  eye(str.data.continuum_dim);
        for igauss=1:size(str.postproc.f_e.N,2)
        %-------------------------------------------------
        % sigma decomposition
        %-------------------------------------------------
        F3D                     =  str.grad.F(:,:,igauss);
        J                       =  str.grad.J(igauss);
        pressure                =  trace(sigma(:,:,igauss))/dim;
        sigma_vol_              =  pressure*I;
        sigma_dev_              =  sigma(:,:,igauss) - sigma_vol_;
        %-------------------------------------------------
        % S decomposition
        %-------------------------------------------------
        invF                    =  inv(F3D);
        S                       =  J*(F3D\(sigma(:,:,igauss)*inv(F3D)'));
        S_vol_                  =  J*pressure*invF*invF';
        S_dev_                  =  S - S_vol_;
        %-------------------------------------------------
        % Storage.
        %-------------------------------------------------
        sigma_vol(:,:,igauss)   =  sigma_vol_;
        sigma_dev(:,:,igauss)   =  sigma_dev_;
        p(igauss,1)             =  pressure;
        S_vol(:,:,igauss)       =  S_vol_;
        S_dev(:,:,igauss)       =  S_dev_;
        end
end
