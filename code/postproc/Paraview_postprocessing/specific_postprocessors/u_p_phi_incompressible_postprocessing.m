function [F,H,J,E0,D0,SigmaF,SigmaH,SigmaJ,Kirchhoff,p]  =  u_p_phi_incompressible_postprocessing(str,gauss_level_information,SigmaF,SigmaH,SigmaJ,Kirchhoff,D0,p)

F                                                 =  str.grad.F;
H                                                 =  str.grad.H;
J                                                 =  str.grad.J;
E0                                                =  str.grad.E0;
F3D                                               =  F;
H3D                                               =  H;
switch str.data.dim
    case 2
        F3D(3,3,:)                               =  ones(1,size(F3D,3));
        H3D(3,3,:)                               =  J;
end
for igauss=1:size(str.f_e.N,2)
    DGammaDF                                      =  gauss_level_information.DGammaDF(:,igauss);
    switch str.data.dim
        case 2
            DGammaDF                              =  [DGammaDF(1:2,1); 0;DGammaDF(3:4,1); zeros(3,1);DGammaDF(5,1)];
    end
    DGammaDF                                      =  reshape(DGammaDF,3,3);
    DGammaDF                                      =  DGammaDF';
    Piola1                                        =  DGammaDF;
    DGammaDH                                      =  gauss_level_information.DGammaDH(:,igauss);
    switch str.data.dim 
        case 2
            DGammaDH                              =  [DGammaDH(1:2,1); 0;DGammaDH(3:4,1); zeros(3,1);DGammaDH(5,1)];
    end
    DGammaDH                                      =  reshape(DGammaDH,3,3);
    DGammaDH                                      =  DGammaDH';
    Piola2                                        =  Javier_double_cross_product(DGammaDH,F3D,1,1,3);
    Piola3                                        =  gauss_level_information.DGammaDJ(igauss)*H3D(:,:,igauss);
    Piola                                         =  Piola1 + Piola2 + Piola3;
    Kirchhoff(:,:,igauss)                         =  Piola + gauss_level_information.P(igauss)*H3D(:,:,igauss);
    SigmaF(:,:,igauss)                            =  Piola1;
    SigmaH(:,:,igauss)                            =  DGammaDH;
    SigmaJ(igauss)                                =  gauss_level_information.DGammaDJ(igauss);
    D0(:,igauss)                                  =  -gauss_level_information.DGammaDE0(:,igauss);
    p(igauss)                                     =  gauss_level_information.P(igauss);
end
