%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Newton Raphson algorithm
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function str                      =  NewtonRaphson(str)

while str.NR.accumulated_factor<0.99  
    %----------------------------------------------------------------------
    % Incremental variables for the load increment strategy.               
    %----------------------------------------------------------------------
    str.NR.incr_load              =  str.NR.incr_load + 1; 
    str.NR.accumulated_factor     =  str.NR.accumulated_factor + str.NR.load_factor;
    convergence_warning           =  'desactivated';
    %----------------------------------------------------------------------    
    % Update Dirichlet boundary conditions.   
    %----------------------------------------------------------------------    
    old_solution                  =  str.solution;
    str                           =  UpdateDirichletBoundaryConditions(str);
    %----------------------------------------------------------------------
    % Updated residual for the next load increment taking into account 
    % Dirichlet boundary conditions.       
    %----------------------------------------------------------------------
    str                           =  NewtonRaphsonInitialResidual(str,old_solution);
    %----------------------------------------------------------------------
    % Necessary incremental variables.                               
    %----------------------------------------------------------------------    
    Residual_dimensionless        =  1e8;  
    str.NR.iteration              =  0;
    nonconvergence_criteria       =  Residual_dimensionless>str.NR.tolerance;
    while nonconvergence_criteria    
        tic  
        str.NR.iteration          =  str.NR.iteration + 1;
        %------------------------------------------------------------------
        % solving system of equations                                                                                         
        %------------------------------------------------------------------
        [freedof,fixdof]          =  DeterminationVariableFreeFixedDofs(str);
        str.solution              =  SolveSystemEquations(freedof,fixdof,str.assembly,str.solution);                                                                    
        %------------------------------------------------------------------
        % update the variables of the problem. corrector step.                                
        %------------------------------------------------------------------
        str                       =  FieldsUpdate(str);
        %------------------------------------------------------------------
        % update matrices and force vectors.                                                                                                                                                       
        %------------------------------------------------------------------
        str                       =  TotalAssembly(str);
        str                       =  NewtonRaphsonResidualUpdate(str);
        %------------------------------------------------------------------
        % checking for convergence.                                                                                                                 
        %------------------------------------------------------------------
        [str.NR,Residual_dimensionless,...
            str.assembly]         =  NewtonRaphsonConvergence(str.NR,str.assembly,str.bc);
        %------------------------------------------------------------------
        % screen ouput for current Newton-Raphson iteration.                                 
        %------------------------------------------------------------------
        toc
        NewtonRaphsonIterationPrint(str.NR,Residual_dimensionless);
        %------------------------------------------------------------------
        % Break in case of nonconvergence           
        %------------------------------------------------------------------
        switch convergence_warning
            case 'activated' 
                 str                                  =  old_str;
                break;  
        end  
    end  
    switch convergence_warning
        %------------------------------------------------------------------
        % Non converged results 
        %------------------------------------------------------------------
        case 'activated'      
        %------------------------------------------------------------------
        % Converged results 
        %------------------------------------------------------------------
        case 'desactivated'
            %--------------------------------------------------------------
            % Printing orders for the converged case.    
            %--------------------------------------------------------------
            NewtonRaphsonConvergedSolutionPrint(str.NR);            
            %--------------------------------------------------------------
            % Postprocessing for static case and saving results.                                     
            %--------------------------------------------------------------
            NewtonRaphsonPostprocessing(str);            
            format long e  
    end                   
end                 

 
