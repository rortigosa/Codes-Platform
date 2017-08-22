%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  Initialisation of parameters for the Newton-Raphson algorithm
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function str               =  NRInitialisation(str)  

str.NR.accumulated_factor  =  0;
str.NR.incr_load           =  0;
str.NR.load_factor         =  1/str.NR.n_incr_loads;
