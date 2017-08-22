%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function prints the summary of results after the Newton-Raphson has
% converged for a specific load increment. The printed results are the
% number of iterations needed for convergence and the current acumulated
% load factor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewtonRaphsonConvergedSolutionPrint(NR)

fprintf('\n\n\n--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('summary of results for the current load increment:\n')
fprintf('load increment: %d\n',NR.incr_load)
fprintf('total number of newton-raphson iteration: %d\n',NR.iteration)
fprintf('the load factor is %12.5e\n',NR.load_factor)
fprintf('the total acumulated load factor is %f\n',NR.accumulated_factor)
fprintf('--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n')
fprintf('--------------------------------------------------------------------\n\n\n')
