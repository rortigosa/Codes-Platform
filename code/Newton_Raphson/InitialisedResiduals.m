%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
%  This function initialises residuals
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                   =  InitialisedResiduals(str)

str.assembly.Residual          =  zeros(str.solution.n_dofs,1);
str.assembly.Residual_stored   =  cell(str.NR.n_incr_loads,1);
