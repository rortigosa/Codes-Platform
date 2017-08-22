%------------------------------------------------------
%------------------------------------------------------
% Calculate the energy density function for the final
% formulation in isotropy.
%------------------------------------------------------
%------------------------------------------------------

function [str]                   =  stored_energy_density_isotropic_paper_form(str) 

str.Phi                          =  zeros(size(str.quadrature.Chi,1),1);
str.Phi_mec_iso                  =  str.Phi;
str.Phi_mec_aniso                =  str.Phi;
str.Phi_diel_iso                 =  str.Phi;
str.Phi_diel_aniso               =  str.Phi;
str.Phi_piezo                    =  str.Phi;
m_id                             =  str.properties.material_identifier;
%------------------------------------------------------
% Properties.
%------------------------------------------------------
im1                              =  str.properties.im1(m_id);
im2                              =  str.properties.im2(m_id);
im3                              =  str.properties.im3(m_id);
lambda_iso                       =  str.properties.lambda_iso(:,m_id);
vacuum_permittivity              =  8.854187817620e-12;
c1                               =  max(max(str.properties.im1));
c2                               =  max(max(str.properties.im2));
c3                               =  max(max(str.properties.im3));
c12                              =  max(c1,c2);
c13                              =  max(c1,c3);
cf                               =  max(c12,c13);
c                                =  1/(cf*sqrt(vacuum_permittivity));
for iloop1=1:size(str.quadrature.Chi,1)
    %------------------------------------------------------
    % Recover kinematic variables.
    %------------------------------------------------------
    switch str.data.dim
          case {2,3}
                F                =  str.postproc.F(:,:,iloop1);
           case {12,13}
                F                =  str.postproc.F(:,:,iloop1);
    end
    C                            =  F'*F;
    E0                           =  str.postproc.E0(:,iloop1);
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

    PhiP_diel                    =  PhiP - PhiP_mec_iso;
    %****************************************** 
    % 4. Total contribution. S
    %******************************************
    str.Phi(iloop1,1)            =  PhiP;                          
    str.Phi_mec_iso(iloop1,1)    =  PhiP_mec_iso;
    str.Phi_diel_iso(iloop1,1)   =  PhiP_diel;
end 
str.Phi_mec                      =  str.Phi_mec_iso  + str.Phi_mec_aniso;
str.Phi_diel                     =  str.Phi_diel_iso + str.Phi_diel_aniso;
str.Phi_elec                     =  str.Phi_diel + str.Phi_piezo;

