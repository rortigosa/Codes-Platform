%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Data needed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str                                 =  InitialData(str)

%--------------------------------------------------------------------------
%  Formulation type
%--------------------------------------------------------------------------
str.data.formulation                         =  'u';
str.data.formulation                         =  'up';
str.data.formulation                         =  'FHJ';
str.data.formulation                         =  'CGC';
str.data.formulation                         =  'CGCCascade';
str.data.formulation                         =  'electro_mechanics';
str.data.formulation                         =  'electro_mechanics_Helmholtz';
%str.data.formulation                        =  'electro_mechanics_incompressible';
%str.data.formulation                        =  'electro_mechanics_Helmholtz_incompressible';
%str.data.formulation                        =  'electro_mechanics_mixed_incompressible';
%str.data.formulation                        =  'u_phi_D0_d_Sigmad';
%str.data.formulation                        =  'electro_BEM_FEM';  %  no mechanics!!
%str.data.formulation                        =  'electro';  %  no mechanics!!
%str.data.formulation                        =  'electro_mechanics_BEM_FEM';
%str.data.formulation                        =  'electro_mechanics_Helmholtz_BEM_FEM';
%str.data.formulation                        =  'electro_mechanics_incompressible_BEM_FEM';
%str.data.formulation                        =  'electro_mechanics_Helmholtz_incompressible_BEM_FEM';
%str.data.formulation                        =  'electro_mechanics_mixed_incompressible_BEM_FEM';
%--------------------------------------------------------------------------
%  Static or dynamic analysis
%--------------------------------------------------------------------------
str.data.analysis                            =  'static';
%str.data.analysis                           =  'dynamic';
%--------------------------------------------------------------------------
%  Finite element
%--------------------------------------------------------------------------
str.fem.shape                                =  1;
str.fem.degree.x                             =  2;  %  Finite Element interpolation order for displacements
str.fem.degree.phi                           =  2;  %  Finite Element interpolation order for displacements
str.fem.degree.pressure                      =  1;  % degree of interpolation for the pressure
str.fem.shape                                =  1;  %  0 is triangular (or tetrahedral) and 1 is quadrilateral (or cubic)
str.fem.degree.postprocessing                =  2;
str.fem.degree.F                             =  1;
str.fem.degree.H                             =  1;
str.fem.degree.J                             =  0;
%str.fem.degree.C                            =  1;
%str.fem.degree.G                            =  1;
%str.fem.degree.c                            =  0;
str.fem.degree.D0                            =  1;
str.fem.degree.d                             =  1;
str.fem.degree.F_continuity                  =  'continuous';
str.fem.degree.F_continuity                  =  'discontinuous';
str.fem.degree.H_continuity                  =  'continuous';
str.fem.degree.H_continuity                  =  'discontinuous';
str.fem.degree.J_continuity                  =  'continuous';
str.fem.degree.J_continuity                  =  'discontinuous';
%str.fem.degree.C_continuity                 =  'continuous';
%str.fem.degree.C_continuity                 =  'discontinuous';
%str.fem.degree.G_continuity                 =  'continuous';
%str.fem.degree.G_continuity                 =  'discontinuous';
%str.fem.degree.c_continuity                 =  'continuous';
%str.fem.degree.c_continuity                 =  'discontinuous';
str.fem.degree.D0_continuity                 =  'continuous';
str.fem.degree.D0_continuity                 =  'discontinuous';
str.fem.degree.d_continuity                  =  'continuous';
str.fem.degree.d_continuity                  =  'discontinuous';
switch str.data.formulation
    case {'electro_BEM_FEM','electro_incompressible_BEM_FEM','electro_mixed_incompressible_BEM_FEM'}
         str.fem.degree.q                    =  1;
         str.fem.surface.BEM_FEM.continuity  =  1;
end
%--------------------------------------------------------------------------
%  Gauss integration 
%--------------------------------------------------------------------------
str.quadrature.degree                         =  max(str.fem.degree.x,str.fem.degree.phi);
%--------------------------------------------------------------------------
%  Newton Raphson parameters
%--------------------------------------------------------------------------
str.NR.n_incr_loads                           =  20;
str.NR.tolerance                              =  1e-6;
str.NR.convergence_plotting                   =  1;
%--------------------------------------------------------------------------
%  Time integrator
%--------------------------------------------------------------------------
str.time_integrator.mass_matrix_flag          =  1;
str.time_integrator.type                      =  'Newmark_beta'; 
str.time_integrator.beta                      =  1/4;
str.time_integrator.gamma                     =  1/2;
str.time_integrator.total_time                =  1e-2;
str.time_integrator.n_time_steps              =  100;   % This gives Dt=1e-4 (cfl=1)
str.time_integrator.Deltat                    =  1e-4;
str.time_integrator.time_iteration            =  1;
str.time_integrator.total_time                =  str.time_integrator.total_time*1e2;
str.time_integrator.stopping                  =  1;

str.time_integrator.alpha_viscous             =  0;
str.time_integrator.beta_viscous              =  0;
str.time_integrator.time                      =  0;

%--------------------------------------------------------------------------
%  Contact information
%--------------------------------------------------------------------------
str.contact.lagrange_multiplier               =  0;
str.fem.degree.contact                        =  1;
str.contact.method                            =  'node';



