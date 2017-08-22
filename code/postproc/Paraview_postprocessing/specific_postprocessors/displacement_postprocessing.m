function [F,H,J,SigmaF,SigmaH,SigmaJ,Kirchhoff]   =  displacement_postprocessing(SigmaF,SigmaH,SigmaJ,Kirchhoff,gauss_level_information)

F                                                 =  str.grad.F;
H                                                 =  str.grad.H;
J                                                 =  str.grad.J;
F3D                                               =  F;
H3D                                               =  H;
switch str.data.dim
    case 2
        F3D(3,3,:)                               =  ones(1,size(F3D,3));
        H3D(3,3,:)                               =  J;
end
for igauss=1:size(str.f_e.N,2)
    DWDF                                          =  gauss_level_information.DWDF(:,igauss);
    switch str.data.dim
        case 2
            DWDF                                  =  [DWDF(1:2,1); 0;DWDF(3:4,1); zeros(3,1);DWDF(5,1)];
    end
    DWDF                                          =  reshape(DWDF,3,3);
    DWDF                                          =  DWDF';
    Piola1                                        =  DWDF;
    DWDH                                          =  gauss_level_information.DWDH(:,igauss);
    switch str.data.dim
        case 2
            DWDH                                  =  [DWDH(1:2,1); 0;DWDH(3:4,1); zeros(3,1);DWDH(5,1)];
    end
    DWDH                                          =  reshape(DWDH,3,3);
    DWDH                                          =  DWDH';
    Piola2                                        =  Javier_double_cross_product(DWDH,F3D,1,1,3);
    Piola3                                        =  gauss_level_information.DWDJ(igauss)*H3D(:,:,igauss);
    Piola                                         =  Piola1 + Piola2 + Piola3;
    Kirchhoff(:,:,igauss)                         =  Piola;
    SigmaF(:,:,igauss)                            =  Piola1;
    SigmaH(:,:,igauss)                            =  DWDH;
    SigmaJ(igauss)                                =  gauss_level_information.DWDJ(igauss);
end
