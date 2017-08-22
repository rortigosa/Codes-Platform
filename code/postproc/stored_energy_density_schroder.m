%------------------------------------------------------
%------------------------------------------------------
% Calculate the Second-Piola Kirchhoff stress tensor.
%------------------------------------------------------
%------------------------------------------------------

function [str]                      =  stored_energy_density_schroder(str) 

switch str.data.dim
    case {2,12}
        dim                         =  2;
    case {3,13}
        dim                         =  3;
end
str.Phi                             =  zeros(size(str.quadrature.Chi,1),1);
str.Phi_mec_iso                     =  str.Phi;
str.Phi_mec_aniso                   =  str.Phi;
str.Phi_diel_iso                    =  str.Phi;
str.Phi_diel_aniso                  =  str.Phi;
str.Phi_diel                        =  str.Phi;
str.Phi_piezo                       =  str.Phi;
m_id                                =  str.properties.material_identifier;
%--------------------------------------------------------------------------
% Constant variables needed
%--------------------------------------------------------------------------
N                                   =  str.cons.N(:,m_id);
I                                   =  str.cons.I;
%------------------------------------------------------
% Properties.
%------------------------------------------------------
lambda                              =  str.properties.lambda(m_id);
mu                                  =  str.properties.mu(m_id);
alpha1                              =  str.properties.alpha1(m_id);
alpha2                              =  str.properties.alpha2(m_id);
alpha3                              =  str.properties.alpha3(m_id);
beta1                               =  str.properties.beta1(m_id);
beta2                               =  str.properties.beta2(m_id);
beta3                               =  str.properties.beta3(m_id);
gamma1                              =  str.properties.gamma1(m_id);
gamma2                              =  str.properties.gamma2(m_id);
for iloop1=1:size(str.quadrature.Chi,1)
    %------------------------------------------------------
    % Recover kinematic variables.
    %------------------------------------------------------
    switch str.data.dim
          case {2,3}
                F                   =  str.postproc.F(:,:,iloop1);
           case {12,13}
                F                   =  str.postproc.F(:,:,iloop1);
    end
    C                               =  F'*F;
    E0                              =  str.postproc.E0(:,iloop1);
    %------------------------------------------------------
    % Create necessary variables.
    %------------------------------------------------------   
    switch str.data.dim
        case 2 
          C3D(1:2,1:2)              =  C;
          C3D(3,3)                  =  1;
          E_Lag3D                   =  0.5*(C3D-eye(3));
        case 3
          C3D                       =  C;
          E_Lag3D                   =  0.5*(C3D-I);
    end
    E_Lag                           =  1/2*(C-I);
    %---------------------------------------
    % Invariants.
    %---------------------------------------
    I1                              =  trace(E_Lag3D);
    I2                              =  trace(E_Lag3D*E_Lag3D);
    I4                              =  N'*E_Lag*N;
    I5                              =  N'*E_Lag*E_Lag*N;
    J1                              =  E0'*E0;
    J2                              =  E0'*N; 
    K1                              =  E0'*E_Lag*N;
    %---------------------------------------
    % Tensors.
    %--------------------------------------- 
    %******************************************
    % 4. Total contribution. S
    %******************************************
    PhiP_mec_iso                    =  0.5*lambda*I1^2 + mu*I2;
    PhiP_mec_aniso                  =  alpha1*I5 + alpha2*I4^2 + alpha3*I1*I4;
    PhiP_diel_iso                   =  gamma1*J1;
    PhiP_diel_aniso                 =  gamma2*J2^2;
    PhiP_diel                       =  PhiP_diel_iso + PhiP_diel_aniso;
    PhiP_piezo                      =  beta1*I1*J2 + beta2*J2*I4 + beta3*K1;
    PhiP_mec                        =  PhiP_mec_iso + PhiP_mec_aniso;
    str.Phi(iloop1,1)               =  PhiP_mec + PhiP_piezo;
    str.Phi_mec_iso(iloop1,1)       =  PhiP_mec_iso;
    str.Phi_mec_aniso(iloop1,1)     =  PhiP_mec_aniso;
    str.Phi_piezo(iloop1,1)         =  PhiP_piezo;                     
    str.Phi_diel(iloop1,1)          =  PhiP_diel;
    str.Phi_diel_iso(iloop1,1)      =  PhiP_diel_iso;
    str.Phi_diel_aniso(iloop1,1)    =  PhiP_diel_aniso;
end
str.Phi_mec                         =  str.Phi_mec_iso + str.Phi_mec_aniso;                             
str.Phi_elec                        =  str.Phi_mec + str.Phi_piezo + str.Phi_diel;
