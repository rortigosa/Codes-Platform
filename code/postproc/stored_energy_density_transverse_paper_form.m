%------------------------------------------------------
%------------------------------------------------------
% Calculate the energy density function for the final 
% formulation for transverse isotropy and rhombic
% and orthotropy.
%------------------------------------------------------
%------------------------------------------------------

function [str]                   =  stored_energy_density_transverse_paper_form(str) 

str.Phi                          =  zeros(size(str.quadrature.Chi,1),1);
str.Phi_mec_iso                  =  str.Phi;
str.Phi_mec_aniso                =  str.Phi;
str.Phi_diel_iso                 =  str.Phi;
str.Phi_diel_aniso               =  str.Phi;
str.Phi_piezo                    =  str.Phi;
m_id                             =  str.properties.material_identifier;
%--------------------------------------------------------------------------
% Constant variables needed
%--------------------------------------------------------------------------
N                                =  str.cons.N(:,m_id);
diel_G                           =  str.anisotropy.diel_G(:,:,:,m_id);
piez_G                           =  str.anisotropy.piez_G(:,:,:,m_id);
elast_G_3D                       =  str.anisotropy.elastic_G_3D(:,:,:,m_id);
%------------------------------------------------------
% Properties.
%------------------------------------------------------
im1                              =  str.properties.im1(m_id);
im2                              =  str.properties.im2(m_id);
im3                              =  str.properties.im3(m_id);
mu                               =  str.properties.mu(:,m_id);
lambda                           =  str.properties.lambda(:,m_id);
lambda_iso                       =  str.properties.lambda_iso(:,m_id);
omega                            =  str.properties.omega(:,m_id);
alpha                            =  str.properties.alpha(:,m_id);
beta                             =  str.properties.beta(:,m_id);
gamma                            =  str.properties.gamma(:,m_id);
vacuum_permittivity              =  8.854187817620e-12;
c                                =  1/(max(max(mu))*sqrt(vacuum_permittivity));

for iloop1=1:size(str.quadrature.Chi,1) 
    %------------------------------------------------------
    % Recover kinematic variables.
    %------------------------------------------------------
    switch str.data.dim
          case {2,3}
                F                =  str.grad.F(:,:,iloop1);
           case {12,13}
                F                =  str.grad.F(:,:,iloop1);
    end
    C                            =  F'*F;
    E0                           =  str.grad.E0(:,iloop1);
    E0                           =  c*E0;  % Dimensionless electric field.
    %------------------------------------------------------
    % Create necessary variables. 
    %------------------------------------------------------   
    switch str.data.dim
        case 2 
          C3D(1:2,1:2)           =  C;
          C3D(3,3)               =  1;
          invC                   =  inv(C);
          invC3D(1:2,1:2)        =  invC;
          invC3D(3,3)            =  1;
        case 3
          C3D                    =  C;
          invC                   =  inv(C);
          invC3D                 =  invC;
    end
    %---------------------------------------
    % Invariants.
    %---------------------------------------
    Ic                           =  trace(C3D);
    IIIc                         =  det(C3D);
    CofC                         =  IIIc*invC;
    CofC3D                       =  IIIc*invC3D;
    IIc                          =  trace(CofC3D);
    %******************************************
    % 3. Individual contributions.
    %******************************************
    %------------------------------------------
    %------------------------------------------
    % 3.1 Isotropic component.
    %------------------------------------------
    %------------------------------------------
    PhiP                         =  im1*(Ic - 3) + im2*(IIc - 3)  +  im3*(IIIc - 1) + ...
                                             (im1 + 2*im2 + im3)*(1/IIIc - 1);
    PhiP_mec_iso                 =  PhiP;
    %------------------------------------------
    %------------------------------------------
    % 3.2 Anisotropic component of the mechanical potential.
    %------------------------------------------
    %------------------------------------------    
    for iinvr=1:size(mu,1)
        G_3D                     =  elast_G_3D(:,:,iinvr);
        J4                       =  trace(C3D*G_3D);
        J5                       =  IIIc*trace(invC3D*G_3D);  
        deno1                    =  trace(G_3D)^alpha(iinvr)*(alpha(iinvr) + 1);
        deno2                    =  trace(G_3D)^beta(iinvr)*(beta(iinvr) + 1);
        deno3                    =  gamma(iinvr);              
        PhiP                     =  PhiP + mu(iinvr)*((J4^(alpha(iinvr) + 1) - trace(G_3D)^(alpha(iinvr) + 1))/deno1 + ...
                                                  (J5^(beta(iinvr) + 1) - trace(G_3D)^(beta(iinvr) + 1))/deno2 + ...
                                                  trace(G_3D)*(IIIc^(-gamma(iinvr)) - 1)/deno3);        
    end    
    PhiP_mec_aniso               =  PhiP - PhiP_mec_iso;
    %------------------------------------------
    %------------------------------------------
    % 3.3 Dielectric component. Isotropic part.
    %------------------------------------------
    %------------------------------------------    
    %-------------------------
    % First type.
    %-------------------------
    PhiP                         =  PhiP + lambda_iso(1)*(E0'*(C*E0)     +  (E0'*E0)^2 + Ic^2 + 6/IIIc - 15);
    %-------------------------
    % Second type.
    %-------------------------
    PhiP                         =  PhiP + lambda_iso(2)*(E0'*(CofC*E0)  +  (E0'*E0)^2 + IIc^2 + 12/IIIc - 21);
    %-------------------------
    % Third type.
    %-------------------------
    PhiP                         =  PhiP + lambda_iso(3)*(E0'*E0*Ic    +  (E0'*E0)^2 + Ic^2 + 6/IIIc - 15);
    %-------------------------
    % Forth type.
    %-------------------------
    PhiP                         =  PhiP + lambda_iso(4)*(E0'*E0*IIc   +  (E0'*E0)^2 + IIc^2 + 12/IIIc - 21);
    %-------------------------
    % Fifth type.
    %-------------------------
    PhiP                         =  PhiP + lambda_iso(5)*(E0'*E0*IIIc  +  (E0'*E0)^2 + IIIc^2 + 2/IIIc - 3);

    PhiP_diel_iso                =  PhiP - (PhiP_mec_iso + PhiP_mec_aniso);
    %------------------------------------------
    %------------------------------------------
    % 3.4 Dielectric component. Anisotropic part.
    %------------------------------------------
    %------------------------------------------    
    %-------------------------
    % First type.
    %-------------------------
    U                            =  diel_G(:,:,1)*E0;
    PhiP                         =  PhiP + lambda(1)*(U'*(C*U)     +  (U'*U)^2 + Ic^2 + 6/IIIc - 15);    
    %-------------------------
    % Second type.
    %-------------------------
    U                            =  diel_G(:,:,2)*E0;
    PhiP                         =  PhiP + lambda(2)*(U'*(CofC*U)  +  (U'*U)^2 + IIc^2 + 12/IIIc - 21);
    %-------------------------
    % Third type.
    %-------------------------
    U                            =  diel_G(:,:,3)*E0;
    PhiP                         =  PhiP + lambda(3)*(U'*U*Ic    +  (U'*U)^2 + Ic^2 + 6/IIIc - 15);
    %-------------------------
    % Forth type.
    %-------------------------
    U                            =  diel_G(:,:,4)*E0;
    PhiP                         =  PhiP + lambda(4)*(U'*U*IIc   +  (U'*U)^2 + IIc^2 + 12/IIIc - 21);
    %-------------------------
    % Fifth type.
    %-------------------------
    U                            =  diel_G(:,:,5)*E0;
    PhiP                         =  PhiP + lambda(5)*(U'*U*IIIc  +  (U'*U)^2 + IIIc^2 + 2/IIIc - 3);
 
    PhiP_diel_aniso              =  PhiP - (PhiP_mec_iso + PhiP_mec_aniso + PhiP_diel_iso);  
    %------------------------------------------
    %------------------------------------------
    % 3.5 Piezoelectric component.
    %------------------------------------------    
    %------------------------------------------
    %-------------------------
    % First type.
    %-------------------------
    G                            =  piez_G(:,:,1);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(1)*(U'*(C*U)     +  (U'*U)^2 + Ic^2  + ...
                                                   (3*G*E0 - N)'*(3*G*E0 - N) + N'*(CofC*N) + 7/IIIc +...
                                                   (E0'*E0)^2 - 20); 
    %-------------------------
    % Second type.
    %-------------------------
    G                            =  piez_G(:,:,2);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(2)*(U'*(CofC*U)     +  (U'*U)^2 + IIc^2  + ...
                                                   (3*G*E0 - N)'*(3*G*E0 - N) + N'*(C*N) + 13/IIIc +...
                                                   (E0'*E0)^2 - 26); 
    %-------------------------
    % Third type.
    %-------------------------
    G                            =  piez_G(:,:,3);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(3)*(U'*U*Ic     +  (U'*U)^2 + Ic^2  + ...
                                                (5*G*E0 - N)'*(5*G*E0 - N) + 7/IIIc +...
                                                (E0'*E0)^2 - 21); 
    %-------------------------
    % Forth type. 
    %-------------------------
    G                            =  piez_G(:,:,4);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(4)*(U'*U*IIc     +  (U'*U)^2 + IIc^2  + ...
                                                   (5*G*E0 - N)'*(5*G*E0 - N) + 14/IIIc +...
                                                   (E0'*E0)^2 - 28); 
    %-------------------------
    % Fifth type.
    %-------------------------
    G                            =  piez_G(:,:,5);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(5)*(U'*U*IIIc     +  (U'*U)^2 + IIIc^2  + ...
                                                   (3*G*E0 - N)'*(3*G*E0 - N) + 3/IIIc +...
                                                   (E0'*E0)^2 - 7); 
    %-------------------------
    % Sixth type.
    %-------------------------
    G                            =  piez_G(:,:,6);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(6)*(U'*(C*U)     +  (U'*U)^2 + Ic^2  + ...
                                                   (3*G*E0 + N)'*(3*G*E0 + N) + N'*(CofC*N) + 7/IIIc +...
                                                   (E0'*E0)^2 - 20); 
    %-------------------------
    % Seventh type.
    %-------------------------
    G                            =  piez_G(:,:,7);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(7)*(U'*(CofC*U)     +  (U'*U)^2 + IIc^2  + ...
                                                   (3*G*E0 + N)'*(3*G*E0 + N) + N'*(C*N) + 13/IIIc +...
                                                   (E0'*E0)^2 - 26); 
    %-------------------------
    % Eigth type.
    %-------------------------
    G                            =  piez_G(:,:,8);     
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(8)*(U'*U*Ic     +  (U'*U)^2 + Ic^2  + ...
                                                   (5*G*E0 + N)'*(5*G*E0 + N) + 7/IIIc +...
                                                   (E0'*E0)^2 - 21); 
    %-------------------------
    % Nineth type.  
    %-------------------------
    G                            =  piez_G(:,:,9);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(9)*(U'*U*IIc     +  (U'*U)^2 + IIc^2  + ...
                                                   (5*G*E0 + N)'*(5*G*E0 + N) + 14/IIIc +...
                                                   (E0'*E0)^2 - 28); 
    %-------------------------
    % Tenth type.
    %-------------------------
    G                            =  piez_G(:,:,10);
    U                            =  G*E0 + N;
    PhiP                         =  PhiP + omega(10)*(U'*U*IIIc     +  (U'*U)^2 + IIIc^2  + ...
                                                   (3*G*E0 + N)'*(3*G*E0 + N) + 3/IIIc +...
                                                   (E0'*E0)^2 - 7); 

    PhiP_piezo                   =  PhiP - (PhiP_mec_iso + PhiP_mec_aniso + PhiP_diel_iso + PhiP_diel_aniso);  
    %****************************************** 
    % 4. Total contribution. S
    %******************************************
    str.Phi(iloop1,1)            =  PhiP;                          
    str.Phi_mec_iso(iloop1,1)    =  PhiP_mec_iso;
    str.Phi_mec_aniso(iloop1,1)  =  PhiP_mec_aniso;
    str.Phi_diel_iso(iloop1,1)   =  PhiP_diel_iso;
    str.Phi_diel_aniso(iloop1,1) =  PhiP_diel_aniso;
    str.Phi_piezo(iloop1,1)      =  PhiP_piezo;    
end 
str.Phi_mec                      =  str.Phi_mec_iso  + str.Phi_mec_aniso;
str.Phi_diel                     =  str.Phi_diel_iso + str.Phi_diel_aniso;
str.Phi_elec                     =  str.Phi_diel + str.Phi_piezo;

