%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: flustr2plt_v2(flu,str,filename,opt)
% This routine exports pairs of flu and str to an ASCII plt file
% Usage: flustr2plt_v2(flu,str,filename,opt)
% flu and str are the structures that are going to be exported to the file
% given by "filename". opt selects which variables are written to the file
% Possible options
% For exporting STR
% 'strx': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, ltens]
%
% 'strx0': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the UNdeformed coordinates. 
%       Variables exported: [x0,y0, vx,vy, stress, F, ltens]
%
% 'stradvx': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, FX]
%       where FX=F from eulerian material position
%
% 'strx_xmat': exports str (2D/3D, ISPM/IFEM). 
%       Spatial variables are the deformed coordinates. 
%       Variables exported: [x,y, vx,vy, stress, F, ltens, Xmat, Ymat]
%
% For exporting FLU
% 'flu_avgvx_p_avgf': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Variables exported: [xmed,ymed, avgvx, avgvy, p, avgfx, avgfy]
%
% 'flu_avgvx_p_avgf_vs': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Spatial positions are shared with first zone (concatenate ONLY)
%       Variables exported: [avgvx, avgvy, p, avgfx, avgfy]
%
% 'flu_mz_vx_p_fx': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell EDGES
%       Variables are exported in multiple zones without averaging.
%       Variables exported: {[x,y,p] [x,y,vx,fx],[x,y,vy,fy]]
%
% 'flu_mz_advx': exports flu (2D/3D, MAC). 
%       Spatial variables are the eulerian coordinates of MAC cell centres
%       Edge velocities and forces are averaged to compute cell centre
%       Variables exported: {[xmed,ymed, avgvx, avgvy, p, avgfx, avgfy]
%                             [x,y,xmat],[x,y,ymat]}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [str,flup,strp]                                                   = flustr2plt_v2(flu,str,filename,opt,append,interpolate)


if nargin < 5
    append                                                                 = 0;
end
if nargin < 6
    interpolate                                                            = 0; % by default, we don't interpolate to nodes.
end
header                                                                     = 'Example: mixed cell-centred and point data on 2D/3D MAC/FE mesh in BLOCK format';
str_export                                                                 = 0;
flu_export                                                                 = 0;
fem_export                                                                 = 0;

if isempty(flu)
    switch str.data.dim
        case 2
            flu.control.z                                                  = 0;
        case {3,13,23}
            flu.control.z                                                  = 1;
    end
    flu.control.outputfilename                                             = str.control.outputfilename;
    flu.control.t                                                          = str.time;
end

if ~isfield(flu.control,'z')
    flu.control.z                                                          = 0;
end

if isnumeric(opt)
    switch opt
        case {1,11,12,7,10,14}            
            switch opt                
                case {1,7,11}
                    options                                                = 'strx';
                case 12
                    options                                                = 'strx0';                
                case 10
                    options                                                = 'stradvx';
                case 14
                    options                                                = 'strx_xmat';
            end
        case {2,6,13,9}            
            switch opt
                case {2,6}
                    options                                                = 'flu_avgvx_p_avgf';                
                case 13
                    options                                                = 'flu_mz_vx_p_fx';
                case 9
                    options                                                = 'flu_mz_advx';
            end
    end
else
    options = opt;
end

name_p = {'p'};
name_residual = {'normR'};

if flu.control.z
    name_pos                                                               = {'X','Y','Z'};
    name_xmat                                                              = {'Xmat','Ymat','Zmat'};
    name_avgvel                                                            = {'avgvx','avgvy','avgvz'};
    name_vel                                                               = {'vx','vy','vz'};
    name_avgf                                                              = {'avgfx','avgfy','avgfz'};
    name_force                                                             = {'fx','fy','fz'};
    %---------------------------
    % Second Piola and contributions.
    %---------------------------
    name_stress                                                            = {'sxx','syy','szz','syx','szx','szy'};
    name_stress_volumetric                                                 = {'s_volumetric_xx','s_volumetric_yy','s_volumetric_zz','s_volumetric_yx','s_volumetric_zx','s_volumetric_zy'};
    name_stress_deviatoric                                                 = {'s_deviatoric_xx','s_deviatoric_yy','s_deviatoric_zz','s_deviatoric_yx','s_deviatoric_zx','s_deviatoric_zy'};
    %---------------------------
    % Cauchy and contributions.
    %---------------------------
    name_cauchy                                                            = {'sigmaxx','sigmayy','sigmazz','sigmayx','sigmazx','sigmazy'};
    name_cauchy_volumetric                                                 = {'sigma_volumetric_xx','sigma_volumetric_yy','sigma_volumetric_zz','sigma_volumetric_yx','sigma_volumetric_zx','sigma_volumetric_zy'};
    name_cauchy_deviatoric                                                 = {'sigma_deviatoric_xx','sigma_deviatoric_yy','sigma_deviatoric_zz','sigma_deviatoric_yx','sigma_deviatoric_zx','sigma_deviatoric_zy'};
    name_cauchy_pressure                                                   = {'hydrostatic pressure'};
    %---------------------------
    % Conjugate stresses.
    %---------------------------        
    name_SigmaF                                                            = {'SigmaF11','SigmaF12','SigmaF13','SigmaF21','SigmaF22','SigmaF23','SigmaF31','SigmaF32','SigmaF33'};
    name_SigmaH                                                            = {'SigmaH11','SigmaH12','SigmaH13','SigmaH21','SigmaH22','SigmaH23','SigmaH31','SigmaH32','SigmaH33'};
    name_SigmaJ                                                            = {'SigmaJ'};
    name_First_Piola                                                       = {'First_Piola11','First_Piola12','First_Piola13','First_Piola21','First_Piola22','First_Piola23','First_Piola31','First_Piola32','First_Piola33'};
    %---------------------------
    % Deformation.
    %---------------------------
    name_F                                                                 = {'f11','f12','f13','f21','f22','f23','f31','f32','f33'};
    name_H                                                                 = {'h11','h12','h13','h21','h22','h23','h31','h32','h33'};
    name_J                                                                 = {'J'};
    name_gradJ                                                             =  {'gradJ_1','gradJ_2','gradJ_3'};
    name_ratioJ                                                            =  {'ratioJ'};
    name_ep                                                                =  {'ep'};
    name_lambda_p                                                          =  {'lambda_p1','lambda_p2','lambda_p3'};
    name_FX                                                                = {'fx11','fx12','fx13','fx21','fx22','fx23','fx31','fx32','fx33'};
    name_L                                                                 = {'L11','L12','L13','L21','L22','L23','L31','L32','L33'};
    name_eigval                                                            = {'L1','L2','L3'};
    name_eigvect                                                           = {'V1x','V1y','V1z','V2x','V2y','V2z','V3x','V3y','V3z'};
    name_displ                                                             = {'Ux','Uy','Uz'};
    name_pos0                                                              = {'X0','Y0','Z0'};
    name_vinc                                                              = {'vinc'};
    name_potential                                                         = {'elect_potential'};
    name_E                                                                 = {'eul_elec_field_x','eul_elec_field_y','eul_elec_field_z'};
    name_E0                                                                = {'lag_elec_field_x','lag_elec_field_y','lag_elec_field_z'};
    %---------------------------
    % Eulerian electric displacement.
    %---------------------------
    name_D                                                                 = {'eul_elec_disp_x','eul_elec_disp_y','eul_elec_disp_z'};
    %---------------------------
    % Lagrangian electric displacement.
    %---------------------------
    name_D0                                                                = {'lag_elec_disp_x','lag_elec_disp_y','lag_elec_disp_z'};
    %---------------------------
    % Nodal loads.
    %---------------------------
    name_F_i                                                               = {'F_int_x','F_int_y','F_int_z'};
    %---------------------------
    %Electric intensity.
    %---------------------------
    name_J_electric                                                        = {'eulerian_electric_density_x','eulerian_electric_density_y','eulerian_electric_density_z'};
    name_J0_electric                                                       = {'lagrangian_electric_density_x','lagrangian_electric_density_y','lagrangian_electric_density_z'};
    eltype                                                                 = 'TETRAHEDRON';
else
    name_pos                                                               = {'X','Y'};
    name_xmat                                                              = {'Xmat','Ymat'};
    name_avgvel                                                            = {'avgvx','avgvy'};
    name_vel                                                               = {'vx','vy'};
    name_avgf                                                              = {'avgfx','avgfy'};
    name_force                                                             = {'fx','fy'};
    %---------------------------
    % Second Piola and contributions.
    %---------------------------
    name_stress                                                            = {'sxx','syy','sxy'};
    name_stress_volumetric                                                 = {'s_volumetric_xx','s_volumetric_yy','s_volumetric_yx'};
    name_stress_deviatoric                                                 = {'s_deviatoric_xx','s_deviatoric_yy','s_deviatoric_yx'};
    %---------------------------
    % Cauchy and contributions.
    %---------------------------
    name_cauchy                                                            = {'sigmaxx','sigmayy','sigmayx'};
    name_cauchy_volumetric                                                 = {'sigma_volumetric_xx','sigma_volumetric_yy','sigma_volumetric_yx'};
    name_cauchy_deviatoric                                                 = {'sigma_deviatoric_xx','sigma_deviatoric_yy','sigma_deviatoric_yx'};
    name_cauchy_pressure                                                   = {'hydrostatic pressure'};
    %---------------------------
    % Conjugate stresses.
    %---------------------------        
    name_SigmaF                                                            = {'SigmaF11','SigmaF12','SigmaF21','SigmaF22'};
    name_SigmaH                                                            = {'SigmaH11','SigmaH12','SigmaH21','SigmaH22'};
    name_SigmaJ                                                            = {'SigmaJ'};
    name_First_Piola                                                       = {'First_Piola11','First_Piola12','First_Piola21','First_Piola22'};    
    %---------------------------
    % Deformation .
    %---------------------------
    name_F                                                                 = {'f11','f12','f21','f22'};
    name_H                                                                 = {'h11','h12','h21','h22'};
    name_J                                                                 = {'J'};
    name_gradJ                                                             =  {'gradJ_1','gradJ_2'};
    name_ratioJ                                                            =  {'ratioJ'};
    name_ep                                                                =  {'ep'};
    name_lambda_p                                                          =  {'lambda_p1','lambda_p2'};
    name_L                                                                 = {'L11','L12','L21','L22'};
    name_FX                                                                = {'fx11','fx12','fx21','fx22'};
    name_eigval                                                            = {'L1','L2'};
    name_eigvect                                                           = {'V1x','V1y','V2x','V2y'};
    name_displ                                                             = {'Ux','Uy'};
    name_pos0                                                              = {'X0','Y0'};    
    name_vinc                                                              = {'vinc'};
    name_potential                                                         = {'elect_potential'};
    name_E                                                                 = {'eul_elec_field_x','eul_elec_field_y'};
    name_E0                                                                = {'lag_elec_field_x','lag_elec_field_y'};
    %---------------------------
    % Eulerian electric displacement.
    %---------------------------
    name_D                                                                 = {'eul_elec_disp_x','eul_elec_disp_y'};
    name_D_diel_iso                                                        = {'eul_elec_di_iso_disp_x','eul_elec_di_iso_disp_y'};
    name_D_diel_aniso                                                      = {'eul_elec_di_aniso_disp_x','eul_elec_di_aniso_disp_y'};
    name_D_diel                                                            = {'eul_elec_di_disp_x','eul_elec_di_disp_y'};
    name_D_piezo                                                           = {'eul_elec_piezo_disp_x','eul_elec_piezo_disp_y'};
    %---------------------------
    % Lagrangian electric displacement.
    %---------------------------
    name_D0                                                                = {'lag_elec_disp_x','lag_elec_disp_y'};
    name_D0_diel_iso                                                       = {'lag_elec_di_iso_disp_x','lag_elec_di_iso_disp_y'};
    name_D0_diel_aniso                                                     = {'lag_elec_di_aniso_disp_x','lag_elec_di_aniso_disp_y'};
    name_D0_diel                                                           = {'lag_elec_di_disp_x','lag_elec_di_disp_y'};
    name_D0_piezo                                                          = {'lag_elec_piezo_disp_x','lag_elec_piezo_disp_y'};
    %---------------------------
    % Nodal forces.
    %---------------------------
    name_F_i                                                               = {'F_int_x','F_int_y'};
    %---------------------------
    % Electric intensity vector.
    %---------------------------
    name_J_electric                                                        = {'eulerian_electric_density_x','eulerian_electric_density_y'};
    name_J0_electric                                                       = {'lagrangian_electric_density_x','lagrangian_electric_density_y'};
    eltype                                                                 = 'TRIANGLE';
end

passive                                                                    = {[]};
var_sharing                                                                = {[]};
con_sharing                                                                = {[]};

numzones                                                                   = 1;
size_var                                                                   = {@extract_pressure};


switch options 
    case 'fem_electromech'
        str                                                                = postprocessing(str); 
        str.kindelem.stress                                                = str.postproc.stress;
        nodes                                                              = @(s) s.Eulerian_x;
        nodes0                                                             = @(s) s.Lagrangian_X;        
        switch str.data.example_comparison
            case 0
                 potential                                                 = @(s) s.phi';
                 displ                                                     = @(s) s.Eulerian_x-s.Lagrangian_X;         
            case 1
                 potential                                                 = @(s) s.postproc.phi';
                 displ                                                     = @(s) s.postproc.displacement;         
        end                
        %----------------------------------
        % Different components of the Second Piola.
        %----------------------------------
        S                                                                  = @(s) s.postproc.stress;
        S_mec_iso                                                          = @(s) s.postproc.stress_mec_iso;
        S_mec_aniso                                                        = @(s) s.postproc.stress_mec_aniso;
        S_mec                                                              = @(s) s.postproc.stress_mec;
        S_diel_iso                                                         = @(s) s.postproc.stress_diel_iso;
        S_diel_aniso                                                       = @(s) s.postproc.stress_diel_aniso;
        S_diel                                                             = @(s) s.postproc.stress_diel;
        S_piezo                                                            = @(s) s.postproc.stress_piezo;
        S_elec                                                             = @(s) s.postproc.stress_elec; 
        S_volumetric                                                       = @(s) s.postproc.stress_volumetric;
        S_deviatoric                                                       = @(s) s.postproc.stress_deviatoric;
        %----------------------------------
        % Stress conjugates.
        %----------------------------------
        SigmaF                                                             = @(s) eT(s.postproc.SigmaF);
        SigmaH                                                             = @(s) eT(s.postproc.SigmaH);
        SigmaJ                                                             = @(s) s.postproc.SigmaJ;
        First_Piola                                                        = @(s) eT(s.postproc.First_Piola);
        %----------------------------------
        % Different components of the Cauchy.
        %----------------------------------
        sigma                                                              = @(s) s.postproc.sigma;
        sigma_mec_iso                                                      = @(s) s.postproc.sigma_mec_iso;
        sigma_mec_aniso                                                    = @(s) s.postproc.sigma_mec_aniso;
        sigma_mec                                                          = @(s) s.postproc.sigma_mec;
        sigma_diel_iso                                                     = @(s) s.postproc.sigma_diel_iso;
        sigma_diel_aniso                                                   = @(s) s.postproc.sigma_diel_aniso;
        sigma_diel                                                         = @(s) s.postproc.sigma_diel;
        sigma_piezo                                                        = @(s) s.postproc.sigma_piezo;
        sigma_elec                                                         = @(s) s.postproc.sigma_elec;
        sigma_volumetric                                                   = @(s) s.postproc.sigma_volumetric;
        sigma_deviatoric                                                   = @(s) s.postproc.sigma_deviatoric;
        sigma_pressure                                                     = @(s) s.postproc.sigma_pressure;
        %----------------------------------
        % Electric field (Eulerian and Lagrangian).
        %---------------------------------- 
        E                                                                  = @(s) s.postproc.E;
        E0                                                                 = @(s) s.postproc.E0; 
        %----------------------------------
        % Eulerian electric displacement. 
        %----------------------------------
        D                                                                  = @(s) s.postproc.D;    
        D_diel_iso                                                         = @(s) s.postproc.D_diel_iso;    
        D_diel_aniso                                                       = @(s) s.postproc.D_diel_aniso;    
        D_diel                                                             = @(s) s.postproc.D_diel;    
        D_piezo                                                            = @(s) s.postproc.D_piezo;    
        %----------------------------------
        % Lagrangian electric displacement.
        %----------------------------------
        D0                                                                 = @(s) s.postproc.D0; 
        D0_diel_iso                                                        = @(s) s.postproc.D0_diel_iso;    
        D0_diel_aniso                                                      = @(s) s.postproc.D0_diel_aniso;    
        D0_diel                                                            = @(s) s.postproc.D0_diel;    
        D0_piezo                                                           = @(s) s.postproc.D0_piezo;    
        %----------------------------------
        % Jacobian.
        %----------------------------------
        J                                                                  = @(s) s.postproc.J;
        ratioJ                                                             = @(s) s.postproc.ratioJ;
        gradJ                                                              = @(s) s.postproc.gradJ;
        %----------------------------------
        % Plastic deformation. 
        %----------------------------------
        ep                                                                 = @(s) s.postproc.hardening_variable;
        lambda_p                                                           = @(s) s.postproc.lambda_p;
        %----------------------------------
        % Electric density vector.
        %----------------------------------
        J_electric                                                         = @(s) s.postproc.J_electric;
        J0_electric                                                        = @(s) s.postproc.J0_electric;
        %---------------------------
        % Stored energy function.
        %---------------------------
        Phi                                                                = @(s) s.postproc.Phi;
        Phi_mec_iso                                                        = @(s) s.postproc.Phi_mec_iso;
        Phi_mec_aniso                                                      = @(s) s.postproc.Phi_mec_aniso;
        Phi_mec                                                            = @(s) s.postproc.Phi_mec;
        Phi_diel_iso                                                       = @(s) s.postproc.Phi_diel_iso;
        Phi_diel_aniso                                                     = @(s) s.postproc.Phi_diel_aniso;
        Phi_diel                                                           = @(s) s.postproc.Phi_diel;
        Phi_piezo                                                          = @(s) s.postproc.Phi_piezo;
        Phi_elec                                                           = @(s) s.postproc.Phi_elec;
        %----------------------------------
        % Force vector.
        %----------------------------------
        switch str.data.dim
            case {2,3}
                dim                                                        = str.data.dim;
                F_i                                                        = @(s) reshape(s.assemb_force.total_force(1:(dim*str.n_nodes)),dim,str.n_nodes);
            case {13,23}
                n_mechanical_nodes                                         = size(str.solid.BEAM_SHELL.discrete.Lagrangian_X,2);
                dim                                                        = 3; 
                F_i                                                        = @(s) zeros(dim,str.n_nodes);
        end
        %----------------------------------
        % Deformation gradient.
        %----------------------------------
        F                                                                  = @(s) eT(s.postproc.F);
        H                                                                  = @(s) eT(s.postproc.H);
        eigenvalsandvects                                                  = @(s) comp_eigenvals_and_eigenvects(s);
        eldata                                                             = @(s) s.connectivity;
        %----------------------------------
        % Data.
        %----------------------------------        
        data                                                               = {{nodes,nodes0,potential,displ,S,S_mec_iso,S_mec_aniso,S_mec,S_diel_iso,S_diel_aniso,S_diel,S_piezo,S_elec,S_volumetric,S_deviatoric,...
                                                                               sigma,sigma_mec_iso,sigma_mec_aniso,sigma_mec,sigma_diel_iso,sigma_diel_aniso,sigma_diel,sigma_piezo,sigma_elec,sigma_volumetric,...
                                                                               sigma_deviatoric,sigma_pressure,SigmaF,SigmaH,SigmaJ,First_Piola,E,E0,D,D_diel_iso,D_diel_aniso,D_diel,D_piezo,D0,D0_diel_iso,D0_diel_aniso,D0_diel,D0_piezo,...
                                                                               F_i,F,H,J,ratioJ,gradJ,ep,lambda_p,J_electric,J0_electric,Phi,Phi_mec_iso,Phi_mec_aniso,Phi_mec,Phi_diel_iso,Phi_diel_aniso,Phi_diel,Phi_piezo,...
                                                                               Phi_elec,eigenvalsandvects, eldata}};
        %----------------------------------
        % Names of the variables associated to data.
        %----------------------------------        
        type                                                               = ones(numel(data{1}),1);
        type(end)                                                          = 0;
        variables                                                          = [name_pos name_pos0 name_potential name_displ name_stress name_stress_mec_iso name_stress_mec_aniso...
                                                                              name_stress_mec name_stress_diel_iso name_stress_diel_aniso name_stress_diel name_stress_piezo name_stress_elec...
                                                                              name_stress_volumetric name_stress_deviatoric name_cauchy name_cauchy_mec_iso name_cauchy_mec_aniso name_cauchy_mec...
                                                                              name_cauchy_diel_iso name_cauchy_diel_aniso name_cauchy_diel name_cauchy_piezo name_cauchy_elec name_cauchy_volumetric...
                                                                              name_cauchy_deviatoric name_cauchy_pressure name_SigmaF name_SigmaH name_SigmaJ name_First_Piola name_E name_E0 name_D name_D_diel_iso name_D_diel_aniso name_D_diel name_D_piezo...
                                                                              name_D0 name_D0_diel_iso name_D0_diel_aniso name_D0_diel name_D0_piezo name_F_i name_F name_H name_J name_ratioJ name_gradJ name_ep name_lambda_p name_J_electric name_J0_electric...
                                                                              name_Phi name_Phi_mec_iso name_Phi_mec_aniso name_Phi_mec name_Phi_diel_iso name_Phi_diel_aniso ...
                                                                              name_Phi_diel name_Phi_piezo name_Phi_elec name_eigval name_eigvect];        
        fem_export                                                          = 1;
        switch str.data.dim
            case 2
                switch str.data.shape
                    case 1 % quadrilateral
                        eltype                                             = 'QUADRILATERAL';
                    case 0 % triangular
                        eltype                                             = 'TRIANGLE';
                end
            case {3,13,23}
                switch str.data.shape
                    case 1 % hexahedral
                        eltype                                             = 'BRICK';
                    case 0 % tetrahedral
                        eltype                                             = 'TETRAHEDRON';
                end
            
        end
        datatype = 'FEBLOCK';

        
    case 'dgfem_u'
        str                                                                = drawcw2Du_p_output(str,str.U);
        nodes                                                              = @(s) s.postproc.nodes; %s.mesh.p;
        displ                                                              = @(s) s.postproc.displ;
        S                                                                  = @(s) s.postproc.stress;
        eldata                                                             = @(s) s.postproc.t'; % s.mesh.t';
        data                                                               = {{nodes,displ,S,eldata}};
        type                                                               = ones(numel(data{1}),1);
        type(end)                                                          = 0;
        variables                                                          = [name_pos name_displ name_stress];        
        fem_export                                                         = 1;
        datatype                                                           = 'FEBLOCK';

    case 'dgfem_u_residual'
        str                                                                = drawcw2Du_p_output(str,str.U);
        nodes                                                              = @(s) s.postproc.nodes; %s.mesh.p;
        displ                                                              = @(s) s.postproc.displ;
        S                                                                  = @(s) s.postproc.stress;
        eldata                                                             = @(s) s.postproc.t'; % s.mesh.t';
        resi                                                               = @(s) s.postproc.residual;
        data                                                               = {{nodes,displ,resi,S,eldata}};
        type                                                               = ones(numel(data{1}),1);
        type(end)                                                          = 0;
        variables                                                          = [name_pos name_displ name_residual name_stress];        
        fem_export                                                         = 1;
        datatype                                                           = 'FEBLOCK';
    
    case 'dgfem_u_residual_only'
        str                                                                = drawcw2Du_p_output(str,str.U);
        nodes                                                              = @(s) s.postproc.nodes; %s.mesh.p;
        displ                                                              = @(s) s.postproc.displ;
        S                                                                  = @(s) s.postproc.stress;
        eldata                                                             = @(s) s.postproc.t'; % s.mesh.t';
        resi                                                               = @(s) s.postproc.residual;
        data                                                               = {{nodes,displ,resi,eldata}};
        type                                                               = ones(numel(data{1}),1);
        type(end)                                                          = 0;
        variables                                                          = [name_pos name_displ name_residual];        
        fem_export                                                         = 1;
        datatype                                                           = 'FEBLOCK';
        
    case {'strx','strx0','stradvx','strx_xmat','strx_cs','strx0_csvs','strx0_eigs','strx_x0_eigs','strx_x0_eigs_vinc','strx_x0_eigs_cs','strx_x0_eigs_csvs','strx0_eigs_cs','strx0_eigs_csvs','strx_eigs','strx_eigs_cs'} % strs (2D/2D ISPM/IFEM w/wo matX)
        nodes                                                              = @(s) s.kindelem.x;
        vel                                                                = @(s) s.kindelem.veloc;       
        S                                                                  = @(s) s.kindelem.stress;
        F                                                                  = @(s) eT(s.kindelem.ftens);
        L                                                                  = @(s) eT(s.kindelem.ltens);
        eldata                                                             = @(s) s.kindelem.lnods';
        xmat                                                               = @(s) s.x;
        eigenvals                                                          = @(s) comp_eigenvals(s);
        eigenvects                                                         = @(s) comp_eigenvects(s);
        eigenvalsandvects                                                  = @(s) comp_eigenvals_and_eigenvects(s);
        str_export                                                         = 1;        
        
        data                                                               = {{nodes,vel,S,F,L,eldata}};
        use_peeked                                                         = ones(numel(data{1}),1);
        use_peeked(1:2)                                                    = 0;
        type                                                               = ones(numel(data{1}),1);
        type(end)                                                          = 0;
                
        variables                                                          = [name_pos name_vel name_stress name_F name_L];
        % All the above should be seen as "defaults"
        
        % We patch them below if necessary
        switch options
            case 'strx'
                nodes                                                      = @(s) s.kindelem.x;
            case 'strx0'                
                nodes                                                      = @(s) s.x;
            case 'stradvx'
                variables                                                  = [name_pos name_vel name_stress name_F name_FX];                                
                L                                                          = @(s) eT(s.kindelem.ftensx);
            case 'strx_xmat'
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_xmat];
                % UNFINISHED. I need to find a way to modify the number of
                % outputs done in strs                
                data                                                       = {{nodes,vel,S,F,L,xmat,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
                type(end)                                                  = 0; % last var is %d
            case 'strx0_eigs'
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) s.x;
                data                                                       = {{nodes,vel,S,F,L,eigenvalsandvects,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
                type(end)                                                  = 0; % last var is %d
            case {'strx_x0_eigs','strx_x0_eigs_csvs','strx_x0_eigs_cs'}
                variables                                                  = [name_pos name_pos0 name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) s.kindelem.x;
                nodes0                                                     = @(s) s.x;
                data                                                       = {{nodes,nodes0,vel,S,F,L,eigenvalsandvects,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:3)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
                type(end)                                                  = 0; % last var is %d
            case 'strx_x0_eigs_vinc'
                variables                                                  = [name_pos name_pos0 name_vel name_stress name_F name_L name_eigval name_eigvect name_vinc];
                nodes                                                      = @(s) s.kindelem.x;
                nodes0                                                     = @(s) s.x;
                vinc                                                       = @(s) s.kindelem.vinc;
                data                                                       = {{nodes,nodes0,vel,S,F,L,eigenvalsandvects,vinc,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:3)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
                type(end)                                                  = 0; % last var is %d

            case 'strx_eigs'
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) s.kindelem.x;
                data                                                       = {{nodes,vel,S,F,L,eigenvalsandvects,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
                type(end)                                                  = 0; % last var is %d    
        end        
        
        switch options
            case {'strx_cs'}
                con_sharing                                                = {'1'};
                data                                                       = {{nodes,vel,S,F,L}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);                
            case {'strx0_csvs'}
                con_sharing                                                = {'1'};
                if flu.control.z
                    var_sharing                                            = {'[1-3]=1'};
                else
                    var_sharing                                            = {'[1-2]=1'};
                end
                data                                                       = {{vel,S,F,L}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1)                                              = 0;
                type                                                       = ones(numel(data{1}),1);                
            case {'strx_x0_eigs_csvs'}
                con_sharing                                                = {'1'};
                if flu.control.z
                    var_sharing                                            = {'[4-6]=1'};
                else
                    var_sharing                                            = {'[3-4]=1'};
                end
                %variables                                                 = [name_pos name_pos0 name_vel name_stress name_F name_L name_eigval name_eigvect];                
                data                                                       = {{nodes,vel,S,F,L,eigenvalsandvects}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);                               
            case {'strx_x0_eigs_cs'}
                con_sharing                                                = {'1'};
                %if flu.control.z
                %    var_sharing                                           = {'[4-6]=1'};
                %else
                %    var_sharing                                           = {'[3-4]=1'};
                %end
                %variables                                                 = [name_pos name_pos0 name_vel name_stress name_F name_L name_eigval name_eigvect];                
                data                                                       = {{nodes,nodes0,vel,S,F,L,eigenvalsandvects,eldata}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:3)                                            = 0;
                type                                                       = ones(numel(data{1}),1);                                               
             case 'strx0_eigs_cs'
                con_sharing                                                = {'1'}; 
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) s.x;
                data                                                       = {{nodes,vel,S,F,L,eigenvalsandvects}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);                
%                 type(end)                                                = 0; % last var is %d 
             case 'strx_eigs_cs'
                con_sharing                                                = {'1'}; 
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) s.kindelem.x;
                data                                                       = {{nodes,vel,S,F,L,eigenvalsandvects}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1:2)                                            = 0;
                type                                                       = ones(numel(data{1}),1);
%                 type(end)                                                = 0; % last var is %d                 
             case 'strx0_eigs_csvs'
                con_sharing                                                = {'1'};
                if flu.control.z
                    var_sharing                                            = {'[1-3]=1'};
                else
                    var_sharing                                            = {'[1-2]=1'};
                end
                variables                                                  = [name_pos name_vel name_stress name_F name_L name_eigval name_eigvect];
                nodes                                                      = @(s) str.x;
                data                                                       = {{vel,S,F,L,eigenvalsandvects}};
                use_peeked                                                 = ones(numel(data{1}),1);
                use_peeked(1)                                              = 0;
                type                                                       = ones(numel(data{1}),1);
%                 type(end)                                                = 0; % last var is %d
           % otherwise  
        end        
        datatype                                                           = 'FEBLOCK';
        
    case {'flu_avgvx_p_avgf','flu_avgvx_p_avgf_vs'} % mac flu 2D/3D       
        variables                                                          = [name_pos name_avgvel name_p name_avgf];        
        nodes                                                              = @(f) ext_fluid_midgrid(f);
        veloc                                                              = @(f) mctmat(avg(ext_fluid_vel(f)));
        p                                                                  = @(f) mctmat({extract_pressure(f)});
        force                                                              = @(f) mctmat(avg(ext_forces(f)));
        flu_export                                                         = 1;
        data                                                               = {{nodes,veloc,p,force}};
        datatype                                                           = 'BLOCK';
        switch  options
            case 'flu_avgvx_p_avgf_vs'
                data                                                       = {{veloc,p,force}};
                if flu.control.z
                    var_sharing                                            = {'[1-3]=1'};
                else
                    var_sharing                                            = {'[1-2]=1'};
                end
        end
        
    case 'flu_mz_vx_p_fx' % multizone mac flu 2D/3D without averaging, including ghosts
        % I don't like this option because tecplot is unable to plot
        % velocity vectors and lots of zone manipulation is required to
        % just show contour levels (the problem is that tecplot assumes a
        % zero value for passive variables in a zone
        variables                                                          = [name_pos name_vel name_p name_force];
        nodes                                                              = @(f) ext_fluid_midgrid(f);
        p                                                                  = @(f) mctmat({extract_pressure(f)});
        nodesX                                                             = @(f) ext_fluid_xgrid(f);
        nodesY                                                             = @(f) ext_fluid_ygrid(f);
        nodesZ                                                             = @(f) ext_fluid_zgrid(f);
        velx                                                               = @(f) mctmat({f.vx});
        vely                                                               = @(f) mctmat({f.vy});
        velz                                                               = @(f) mctmat({f.vz});
        forcex                                                             = @(f) mctmat({f.fx});
        forcey                                                             = @(f) mctmat({f.fy});
        forcez                                                             = @(f) mctmat({f.fz});
        flu_export                                                         = 1;
        datatype                                                           = 'BLOCK';
        
        if flu.control.z
            data                                                           = {{nodes, p}, {nodesX, velx, forcex}, {nodesY, vely, forcey} , {nodesZ, velz, forcez}};
            size_var                                                       = {@extract_pressure, @(f) f.vx, @(f) f.vy, @(f) f.vz};
            numzones                                                        4;
            passive                                                        = {'4-6,8-10','5-7,9-10','4,6-8,10','4-5,7-9'};
        else
            data                                                           = {{nodes, p}, {nodesX, velx, forcex}, {nodesY, vely, forcey}};
            size_var                                                       = {@extract_pressure, @(f) f.vx, @(f) f.vy};
            numzones                                                       = 3;
            passive                                                        = {'3-4,6-7','4-5,7','3,5,6'};
        end
        
    case 'flu_mz_advx' % multizone mac flu 2D/3D with X and Y (material positions as spatial variables)
        variables                                                          = [name_pos name_avgvel name_p name_avgf name_xmat];
        nodes                                                              = @(f) ext_fluid_midgrid(f);
        nodesX                                                             = @(f) ext_fluid_xgrid(f);
        nodesY                                                             = @(f) ext_fluid_ygrid(f);
        nodesZ                                                             = @(f) ext_fluid_zgrid(f);
        Xmat                                                               = @(f) mctmat({f.X});
        Ymat                                                               = @(f) mctmat({f.Y});
        Zmat                                                               = @(f) mctmat({f.Z});
        veloc                                                              = @(f) mctmat(avg(ext_fluid_vel(f)));
        p = @(f) mctmat({extract_pressure(f)});
        force                                                              = @(f) mctmat(avg(ext_forces(f)));
        flu_export                                                         = 1;
                
        datatype = 'BLOCK';
        if flu.control.z
            numzones                                                       = 4;
            data                                                           = {{nodes,veloc,p,force},{nodesX, Xmat},{nodesY, Ymat},{nodesZ, Zmat}};
            size_var                                                       = {@extract_pressure, @(f) f.X, @(f) f.Y, @(f) f.Z};
            passive                                                        = {'11-13', '4-10,12-13','4-11,13','4-12'};
        else
            numzones                                                       = 3;
            data                                                           = {{nodes,veloc,p,force},{nodesX, Xmat},{nodesY, Ymat}};
            size_var                                                       = {@extract_pressure, @(f) f.X, @(f) f.Y};
            passive                                                        = {'8-9', '3-7,9','3-7,8'};
        end
end

loc                                                                        = 'NODAL';
loc                                                                        = 'CELLCENTERED';

if str_export 
switch str.kindelem.eltype
    case {'3dtetra4_nodal','2dtria3nodal_opt','2dtria3nodal_gq','3dtetra4_nodal_gq'}
        ifem                                                               = 0; % I need to finish detecting the type
    otherwise
        ifem                                                               = 1;
end
else
    ifem                                                                   = 0;
end

% Variable location depends on the method. For the ISPM all str-exported
% variables are nodal, whereas for the ifem node positions and velocities
% are nodal, but the rest of the variables are cell-centered.
% For mac flus variables are exported as "nodal", even though strictly
% speaking they are not.
if (ifem && ~flu_export)
    location                                                               = cellfun(@(t) 'CELLCENTERED', variables,'UniformOutput',false);
    mult                                                                   = 2;
    switch options
        case {'strx_x0_eigs','strx_x0_eigs_cs','strx_x0_eigs_csvs'}
            mult                                                           = 3;
    end
    for i = 1:(mult*size(str.x,1))
        location{i}                                                        = 'NODAL';
    end
    switch options
        case 'strx_x0_eigs_vinc'
            location{end}                                                  = 'NODAL';
    end
else
    location                                                               = cellfun(@(t) 'NODAL', variables,'UniformOutput',false);
%     switch options
%         case 'dgfem_u_residual'
%             switch str.d
%                 case 2
%                     location{5}                                          = 'CELLCENTERED';    
%                 case 3
%                     location{7}                                          = 'CELLCENTERED';    
%             end
%     end
end

if isunix()
    dirsep                                                                 = '/';
else
    dirsep                                                                 = '\';
end
a=find(flu.control.outputfilename==dirsep);
if isunix()
    zonenumber                                                             = flu.control.outputfilename((a(end-2)+1):(a(end)-1));
    zonenumber                                                             = strrep(zonenumber,dirsep,'\\');
else
    if ~isempty(a)
        zonenumber                                                         = flu.control.outputfilename(1:(a(end)-1));
        zonenumber                                                         = strrep(zonenumber,dirsep,'\\');
    else
        zonenumber                                                         = 'not automatically assigned';
    end
end



switch str.data.example_comparison
       case 0
            jobfolder                                                      = [str.jobfolder '\results'];
            cd(jobfolder)
       case 1            
            cd(str.destination_jobfolder)
end

%datatype = 'FEPOINT';
% Open file for writing

if append
    f                                                                      = fopen(filename,'a+');
else
    f                                                                      = fopen(filename,'w+');
end

% Write header
writeheader(f,header);

% Write variable header
writevariableheader(f,variables);

if fem_export
     for zon = 1:numzones
        % Write zone header
        writezoneheader_t(f,zonenumber, size(str.postproc.nodes,2),size(str.postproc.t,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 

        % Write variable location
        writevarlocation(f,location);
        
        % Now get the data ready to output
    % STR ISPM
        % In the case of a structure
%         if ispm_str
%             [flup,strp]                                                  = peek_structure(flu,str,15);
%         end
        % write node positions and velocities; always in block format and
        % nodal
        
        % Now get the data ready to output
        for i = 1:numel(data{zon})
            switch type(i)
                case 1
                    formato                                                = '%.16g';
                case 0
                    formato                                                = '%d';
            end
%             switch use_peeked(i)
%                 case 1
%                     s                                                    = strp;
%                 case 0
                    s                                                      = str;
%             end
            writepointdata_faster(f,data{zon}{i}(s),formato);
        end 
        

     end
end

if str_export
     for zon = 1:numzones
             % STR ISPM
        % In the case of a structure
%         if ispm_str
            [flup,strp]                                                    = peek_structure(flu,str,15);
%         end
        % write node positions and velocities; always in block format and
        % nodal
        if interpolate
            if sum(use_peeked) > 0
                str_interp                                                 = resample_str_to_nodes(strp);
            else
                error('Cannot interpolate not peeking the structure!!')
            end
            writezoneheader_t(f,zonenumber, size(str_interp.kindelem.x,2),size(str_interp.kindelem.lnods,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 
        else
            writezoneheader_t(f,zonenumber, size(str.kindelem.x,2),size(str.kindelem.lnods,2),flu.control.t,datatype,eltype,passive{zon},var_sharing{zon},con_sharing{zon}); 
        end
        
        % Write zone header
        
        % Write variable location
        writevarlocation(f,location);
        
        % Now get the data ready to output

        
        % Now get the data ready to output
        for i = 1:numel(data{zon})
            switch type(i)
                case 1
                    formato                                                = '%.16g';
                case 0
                    formato                                                = '%d';
            end
            switch use_peeked(i)
                case 1
                    s                                                      = strp;
                case 0
                    s                                                      = str;
            end
            if interpolate
                s                                                          = str_interp;
            end
            writepointdata_faster(f,data{zon}{i}(s),formato);
        end
        
% %         blockdata = [nodes(str); vel(str)];
% %         % write data
% %         writepointdata_faster(f,blockdata);
% %         % write the other variables; always in block format. nodal or
% %         % cell-centered depending on ispm/ifem
% %         blockdata = [S(strp) ; ...
% %                     F(strp);...
% %                     L(strp)];                
% %         writepointdata_faster(f,blockdata);
% %         % get element data ready
% %         elementdata = (str.kindelem.lnods');
% %         % write element data
% %         writeelementdata_faster(f,elementdata);
% %     % STR IFEM
     end
end

if flu_export
    for zon = 1:numzones
        P                                                                  = size_var{zon}(flu);
        if flu.control.z
            % mac flu 3d
            writezoneheader_ijkt(f,zonenumber, size(P,1),size(P,2), size(P,3),flu.control.t,datatype, passive{zon},var_sharing{zon});
        else
            % mac flu 2d
            writezoneheader_ijt(f,zonenumber, size(P,1),size(P,2),flu.control.t,datatype, passive{zon},var_sharing{zon});
        end 
        % Write variable location
        writevarlocation(f,location);

        % Now get the data ready to output
        for i = 1:numel(data{zon})
            writepointdata_faster(f,data{zon}{i}(flu));
        end
        flup                                                               = flu;
        strp                                                               = str;
    end
end
            
% close file
fclose(f);

end


function writeheader(f,s)
headline                                                                   = ['TITLE = "' s '"'];
fprintf(f,'%s\n',headline);
end
function writevariableheader(f,s)

a                                                                          = cell2mat(cellfun(@(t) ['"' t '", '],s,'UniformOutput',false));
variables                                                                  = ['VARIABLES = ' a(1:(end-2))];
fprintf(f,'%s\n',variables);
end
function writezoneheader(f,zonename,numpoints,numel,type,eltype,varargin)
if nargin > 6
    %passive = 1;  
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive = 0;
    pass                                                                   = '';
end

zones                                                                      = ['ZONE T ="' zonename '", N=' num2str(numpoints) ', E=' num2str(numel) ', F=' type ', ET=' eltype pass ];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ij(f,zonename,numpoints,numel,type,varargin)
if nargin > 6
    %passive                                                               = 1;
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(numpoints) ', J=' num2str(numel) ', F=' type pass];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijt(f,zonename,numpoints,numel,t,type,varargin)
if nargin > 6
    %passive                                                               = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 7
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(numpoints) ', J=' num2str(numel) ', F=' type ', SOLUTIONTIME=' num2str(t,20) pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijk(f,zonename,i,j,k,type, varargin)
if nargin > 6
    %passive = 1;
    pass                                                                   = [', PASSIVEVARLIST=[' varargin{1} ']'];
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 7
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end

zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(i) ', J=' num2str(j) ', K=' num2str(k) ', F=' type pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_ijkt(f,zonename,i,j,k,t,type, varargin)
if nargin > 7
    %passive = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 8
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", I=' num2str(i) ', J=' num2str(j) ', K=' num2str(k) ', F=' type ', SOLUTIONTIME=' num2str(t,20) pass vs];
fprintf(f,'%s\n',zones);
end
function writezoneheader_t(f,zonename,numpoints,numel,t,type,eltype,varargin)
if nargin > 7
    %passive                                                               = 1;
    if ~isempty(varargin{1})
        pass                                                               = [', PASSIVEVARLIST=[' varargin{1} ']'];
    else
        pass                                                               = '';
    end
else
    %passive                                                               = 0;
    pass                                                                   = '';
end
if nargin > 8
    if ~isempty(varargin{2})
        vs                                                                 = [', VARSHARELIST=(' varargin{2} ')'];
    else
        vs                                                                 = '';
    end
else    
    vs                                                                     = '';    
end
if nargin > 9
    if ~isempty(varargin{3})
        cs                                                                 = [', CONNECTIVITYSHAREZONE=' varargin{3}];
    else
        cs                                                                 = '';
    end
else    
    cs                                                                     = '';    
end
zones                                                                      = ['ZONE T ="' zonename '", N=' num2str(numpoints) ', E=' num2str(numel) ', F=' type ', ET=' eltype ', SOLUTIONTIME=' num2str(t,20) pass vs cs];
fprintf(f,'%s\n',zones);
end
function writepointdata(f,elementdata)
% save(filename,'-ascii','-append','-double','pointdata');
c=1;
for i=1:size(elementdata,1)
    for j=1:size(elementdata,2)
        fprintf(f,'%g ',elementdata(i,j));
        c = c+1;
        if mod(c,50)==1 % This is just to avoid extremely long lines and tecplot crashing with them
            fprintf(f,'\n');
        end
    end
    fprintf(f,'\n');
end
end
function writepointdata_fast(filename,elementdata)
save(filename,'-ascii','-append','-double','elementdata');
end
function writeelementdata(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% fprintf(f,'\n');
for i=1:size(elementdata,1)
    for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i,j));
    end
    fprintf(f,'\n');
end
end
function writepointdata_faster(f,elementdata,varargin)
% save(filename,'-ascii','-append','elementdata');
% % fprintf(f,'\n');
if nargin > 2
    formato                                                                = varargin{1};
else
    formato                                                                = '%.16g';
end
stride                                                                     = 100;
elementdata                                                                = elementdata';
 j                                                                         = 0:(stride-1);
 if numel(elementdata) > stride
for i=1:stride:(numel(elementdata)-stride)
%     for j=1:size(elementdata,2)
        fprintf(f,[formato ' '],elementdata(i+j));
%     end
    fprintf(f,'\n');
end
fprintf(f,[formato ' '],elementdata(i+stride:end));
 else
     fprintf(f,[formato ' '],elementdata);
 end
% fprintf(f,'%d ',elementdata'); % this is really fast, but gives long lines that fail with tecplot
fprintf(f,'\n\n');
end
function writeelementdata_fast(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% fprintf(f,'\n');
for i=1:size(elementdata,1)
%     for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i,:));
%     end
    fprintf(f,'\n');
end
end
function writeelementdata_faster(f,elementdata)
% save(filename,'-ascii','-append','elementdata');
% % fprintf(f,'\n');
stride                                                                     = 100;
elementdata                                                                = elementdata';
 j = 0:(stride-1);
 if numel(elementdata) > stride
for i=1:stride:(numel(elementdata)-stride)
%     for j=1:size(elementdata,2)
        fprintf(f,'%d ',elementdata(i+j));
%     end
    fprintf(f,'\n');
end
fprintf(f,'%d ',elementdata(i+stride:end));
 else
     fprintf(f,'%d ',elementdata);
 end
% fprintf(f,'%d ',elementdata'); % this is really fast, but gives long lines that fail with tecplot
fprintf(f,'\n\n');
end
function writedata(filename,elementdata)
save(filename,'-ascii','-append','elementdata');
end
function writevarlocation(f,s)
a                                                                          = cell2mat(cellfun(@(t) [t ', '],s,'UniformOutput',false));
variables                                                                  = ['VARLOCATION = (' a(1:(end-2)) ')'];
fprintf(f,'%s\n',variables);
end
function pointdata                                                         = pointformat(c)
pointdata                                                                  = cell2mat(c);
end
function T                                                                 = eT(tensor)
    if size(tensor,1) == 2
        T                                                                  = [squeeze(tensor(1,1,:))'; ...
                                                                              squeeze(tensor(1,2,:))'; ...
                                                                              squeeze(tensor(2,1,:))'; ...
                                                                              squeeze(tensor(2,2,:))'];
    elseif size(tensor,1) == 3
        T                                                                  = [squeeze(tensor(1,1,:))'; ...
                                                                              squeeze(tensor(1,2,:))'; ...
                                                                              squeeze(tensor(1,3,:))'; ...
                                                                              squeeze(tensor(2,1,:))'; ...
                                                                              squeeze(tensor(2,2,:))'; ...
                                                                              squeeze(tensor(2,3,:))'; ...
                                                                              squeeze(tensor(3,1,:))'; ...
                                                                              squeeze(tensor(3,2,:))'; ...
                                                                              squeeze(tensor(3,3,:))'];
    end
end
function v                                                                 = ext_fluid_vel(f)
    if f.control.z
        vx                                                                 = f.vx(f.index1,f.pindex2,f.pindex3);
        vy                                                                 = f.vy(f.pindex1,f.index2,f.pindex3);
        vz                                                                 = f.vz(f.pindex1,f.pindex2,f.index3);
        v                                                                  = {vx;vy;vz};
    else
        vx                                                                 = f.vx(f.index1,f.pindex2);
        vy                                                                 = f.vy(f.pindex1,f.index2);
        v                                                                  = {vx;vy};
    end
end
function v= ext_forces(f)
    if f.control.z
        vx                                                                 = f.fx(f.index1,f.pindex2,f.pindex3);
        vy                                                                 = f.fy(f.pindex1,f.index2,f.pindex3);
        vz                                                                 = f.fz(f.pindex1,f.pindex2,f.index3);
        v                                                                  = {vx;vy;vz};
    else                                                                   
        vx                                                                 = f.fx(f.index1,f.pindex2);
        vy                                                                 = f.fy(f.pindex1,f.index2);
        v                                                                  = {vx;vy};
    end   
end
function p                                                                 = extract_pressure(f)
     if f.control.z
         p                                                                 = f.p(f.pindex1,f.pindex2,f.pindex3);
     else
         p                                                                 = f.p(f.pindex1,f.pindex2);
     end
     %p={p};
end
function po = mctmat(p)
    switch numel(p)
        case 3
            avx                                                            = p{1}; avy = p{2}; avz = p{3};
            po                                                             = [avx(:) avy(:) avz(:)]';
        case 2
            avx                                                            = p{1}; avy = p{2};
            po                                                             = [avx(:) avy(:)]';
        case 1
            avx                                                            = p{1};
            po                                                             = [avx(:)]';
    end
end
function vo                                                                = avg(v)
    if numel(v)==3
        vx  = v{1}; vy = v{2}; vz = v{3};
        avx                                                                = 0.5*(vx(1:end-1,:,:)+vx(2:end,:,:));
        avy                                                                = 0.5*(vy(:,1:end-1,:)+vy(:,2:end,:));
        avz                                                                = 0.5*(vz(:,:,1:end-1)+vz(:,:,2:end));
        vo                                                                 = {avx;avy;avz};        
    else
        vx                                                                 = v{1}; vy = v{2};% vz = v{3};
        avx                                                                = 0.5*(vx(1:end-1,:)+vx(2:end,:));
        avy                                                                = 0.5*(vy(:,1:end-1)+vy(:,2:end)); 
        vo                                                                 = {avx;avy};
    end    
end
function xo                                                                = ext_fluid_midgrid(f)
    if f.control.z
        % CORRECT THIS IF NECESSARY FOR 3D
        % uncomment if there is any problem
%         xmed=((f.x(1:end-1,1:end-1,1:end-1)+f.x(2:end,1:end-1,1:end-1))/2);
%         ymed=((f.y(1:end-1,1:end-1,1:end-1)+f.y(1:end-1,2:end,1:end-1))/2);
%         zmed=((f.z(1:end-1,1:end-1,1:end-1)+f.z(1:end-1,1:end-1,2:end))/2);
%         xo = [xmed(:) ymed(:) zmed(:)]';
        xo                                                                 = [f.xmed(:) f.ymed(:) f.zmed(:)]';
    else        
        xmed                                                               = ((f.x(1:end-1,1:end-1)+f.x(2:end,1:end-1))/2);
        ymed                                                               = ((f.y(1:end-1,1:end-1)+f.y(1:end-1,2:end))/2);
        xo                                                                 = [xmed(:) ymed(:)]';        
    end
end
function xo                                                                = ext_fluid_xgrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                  (f.xmin(1)-f.deltax*f.nxghostcells):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells),...
                                                                                  (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay),...
                                                                                  (f.xmin(3)-f.deltaz*f.nzghostcells+0.5*f.deltaz):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells-0.5*f.deltaz));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        [xr,yr]                                                            = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay));
        xo                                                                 = [xr(:) yr(:)]';                
    end
end
function xo                                                                = ext_fluid_ygrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells),...
                                                                                    (f.xmin(3)-f.deltaz*f.nzghostcells+0.5*f.deltaz):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells-0.5*f.deltaz));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        [xr,yr]                                                            = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells));
        xo = [xr(:) yr(:)]';                
    end
end
function xo                                                                = ext_fluid_zgrid(f)
    if f.control.z        
        [xr,yr,zr]                                                         = ndgrid(...
                                                                                    (f.xmin(1)-f.deltax*f.nxghostcells+0.5*f.deltax):f.deltax:(f.xmax(1)+f.deltax*f.nxghostcells-0.5*f.deltax),...
                                                                                    (f.xmin(2)-f.deltay*f.nyghostcells+0.5*f.deltay):f.deltay:(f.xmax(2)+f.deltay*f.nyghostcells-0.5*f.deltay),...
                                                                                    (f.xmin(3)-f.deltaz*f.nzghostcells):f.deltaz:(f.xmax(3)+f.deltaz*f.nzghostcells));
        xo                                                                 = [xr(:) yr(:) zr(:)]';                
    else
        xo                                                                 = [];                
    end
end

function e                                                                 = comp_eigenvals_and_eigenvects(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        e                                                                  = zeros(6,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            [V,D]                                                          = eig(A);
            d                                                              = diag(D);
            [e(1:2,i),p]                                                   = sort(d);
            V                                                              = V(:,p);
            e(3:6,i)                                                       = V(:);             
        end
    else
        e                                                                  = zeros(12,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i);... 
                                                                              S(4,i) S(2,i) S(6,i);...
                                                                              S(5,i) S(6,i) S(3,i)];
            [V,D]                                                          = eig(A);
            d                                                              = diag(D);
            [e(1:3,i),p]                                                   = sort(d);
            V                                                              = V(:,p);
            e(4:12,i)                                                      = V(:);
        end
    end
end

function e                                                                 = comp_eigenvals(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        e                                                                  = zeros(2,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    else
        e                                                                  = zeros(3,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i); S(4,i) S(2,i) S(6,i); S(5,i) S(6,i) S(3,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    end
end
function v                                                                 = comp_eigenvects(s)
    S                                                                      = s.kindelem.stress;
    
    if size(S,1) == 3
        v                                                                  = zeros(4,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(3,i); S(3,i) S(2,i)];
            [e,vv]                                                         = sort(eig(A));
        end
    else
        v                                                                  = zeros(3,size(S,2));
        for i = 1:size(S,2)
            A                                                              = [S(1,i) S(4,i) S(5,i); S(4,i) S(2,i) S(6,i); S(5,i) S(6,i) S(3,i)];
            e(:,i)                                                         = sort(eig(A));
        end
    end
end
