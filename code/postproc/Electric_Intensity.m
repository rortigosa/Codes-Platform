function [I,W_dot_elec]                 =  Electric_Intensity(str)

I                                       =  0;
W_dot_elec                              =  0; 
Dt                                      =  str.data.delta_t;
for iedge =1:size(str.surface_elem,2)    
    element_info                        =  str.surface_elem{iedge};
    element                             =  element_info(1);
    nodes                               =  str.connectivity(element,:);
    Normal_vectors                      =  str.flux_calc.N_3D{iedge};
    %-------------------------------------------
    % Material associated. 
    %-------------------------------------------
    str.properties.material_identifier  =  str.properties.material(iedge,1);
    %-------------------------------------------
    % Gauss points and shape functions for the surface considered.
    %-------------------------------------------
    str.quadrature.Chi                  =  [];
    str.f_e.N                           =  [];
    str.f_e.DN_chi                      =  [];
    str.quadrature.Chi                  =  str.flux_calc.quadrature.Chi_3D{iedge};
    str.f_e.N                           =  str.flux_calc.f_e.N_3D{iedge};
    str.f_e.DN_chi                      =  str.flux_calc.f_e.DN_chi_3D{iedge};
    N_shape                             =  str.f_e.N;
%   DN_chi                              =  str.f_e.DN_chi;
    %------------------------ 
    % Element information in 
    % the previous time step.
    %------------------------
    xelem                               =  str.Eulerian_xn1(:,nodes); % Current time step.
    Xelem                               =  str.Lagrangian_X(:,nodes); % Current time step.
    phielem_n1                          =  str.phin1(nodes);          % Current time step.
    %--------------------------------------------------------------------------
    % D0 in previous time step.
    %--------------------------------------------------------------------------
    [str]                               =  gradients(xelem,Xelem,phielem_n1,str);
    Fn1                                 =  str.grad.F;
    Jn1                                 =  str.grad.J; 
    switch str.formulation
        case 'E'
             [D0n1]                     =  electric_displacement(str);
        case 'D_to_E'
             old_str.quadrature         =  str.quadrature;
             str.quadrature.Chi         =  str.flux_calc.quadrature.Chi;
             str.flux_activated         =  1;
             [D0n1,str]                 =  internal_Newton_Raphson(str);                                  
             str.quadrature             =  old_str.quadrature;
             str.flux_activated         =  0;
    end            
    %------------------------
    % Element information in 
    % the current time step.
    %------------------------
    xelem                               =  str.Eulerian_x(:,nodes);   % Current time step.
    Xelem                               =  str.Lagrangian_X(:,nodes); % Current time step.
    phielem                             =  str.phi(nodes);            % Current time step.
    %--------------------------------------------------------------------------
    % D0 in current time step.
    %--------------------------------------------------------------------------
    [str]                               =  gradients(xelem,Xelem,phielem,str);
    Fn                                  =  str.grad.F;
    Jn                                  =  str.grad.J;
    switch str.formulation
        case 'E'
             [D0]                       =  electric_displacement(str);
        case 'D_to_E'
             old_str.quadrature         =  str.quadrature;
             str.quadrature.Chi         =  str.flux_calc.quadrature.Chi;
             str.flux_activated         =  1;
             [D0,str]                   =  internal_Newton_Raphson(str);                                  
             str.quadrature             =  old_str.quadrature;
             str.flux_activated         =  0;
    end                
    %------------------------------------------------------------------
    % Jacobian at every Gauss point.
    %------------------------------------------------------------------
    Jacobian                            =  str.flux_calc.grad.J_3D{iedge};
    %------------------------------------------------------------------
    % Intensity, energy power. main
    %------------------------------------------------------------------    
    for igauss=1:size(str.flux_calc.quadrature.Chi,1)
%       DX_chi                          =  str.grad.DX_chi(:,:,igauss);
%       J_t                             =  DX_chi';
%       J_t                             =  abs(det(J_t));
        J_t                             =  Jacobian(igauss);
        W                               =  str.flux_calc.quadrature.W_v(igauss);
        %---------------------------------------------------------------------
        % Obtain Lagrangian electric displacement.
        %---------------------------------------------------------------------
        D0_                             =  D0(:,igauss);
        D0n1_                           =  D0n1(:,igauss);
        D0Dt                            =  (D0_ - D0n1_)/Dt;
        D_                              =  1/Jn(igauss,1)*Fn(:,:,igauss)*D0(:,igauss);
        Dn1_                            =  1/Jn1(igauss,1)*Fn1(:,:,igauss)*D0n1(:,igauss);
        DDt                             =  (D_ - Dn1_)/Dt;
        phi                             =  phielem'*N_shape(:,igauss);  % Electric potential in the gauss point considered.
        phin1                           =  phielem_n1'*N_shape(:,igauss);  % Electric potential in the gauss point considered.
        DphiDt                          =  (phi - phin1)/Dt;
        Normal_vector                   =  Normal_vectors(:,igauss);
        %-------------------------------------------------------------------
        % Force vector (motion part). 
        %-------------------------------------------------------------------
        I                               =  I  - DDt'*Jn(igauss,1)*(Fn(:,:,igauss)\Normal_vector)*W*J_t;
        W_dot_elec                      =  W_dot_elec  - 0.5*(D0Dt'*Normal_vector*phi + D0_'*Normal_vector*DphiDt)*W*J_t;
    end
end
