function [xopt, fopt, direc, iter, warnflag, funcalls, allvecs] = ...
    powell_conjugate_direction(func, x0, options)
%{
Uses a modification of Powell's method to find the minimum of
a function of N variables with unbounded constaints. 

Outputs:
xopt : array        -Parameter which minimizes 'func'.
fopt : number       -Value of function at minimum: ''fopt = func(xopt)''.
direc : matrix      -Current direction set.
iter : int          -Number of iterations.
funcalls : int      -Number of function calls made.
warnflag : int      -Integer warning flag:
       -1 : User use callback function to terminate
        1 : Maximum number of function evaluations.
        2 : Maximum number of iterations.
        3 : NaN result encountered.
        4 : The result is out of the provided bounds.
allvecs : cell      -List of solutions at each iteration.
%}

arguments
    func function_handle       % Objective function to be optimized.
    x0 (:,1) double            % Initial guess.
    
    % Line search error tolerance
    options.xtol (1,1) double {mustBePositive} = 1e-6
    
    % Relative error in ''func(xopt)'' acceptable for convergence.
    options.ftol (1,1) double {mustBePositive} = 1e-6
    
    % Maximal interations for optimization
    options.maxiter (1,1) double {mustBeInteger, mustBePositive} = length(x0)*10000

    % Maximal times to evaluate objective function
    options.maxfun (1,1) double {mustBeInteger, mustBePositive} = length(x0)*10000

    % If true, print convergence messages.
    options.disp {mustBeNumericOrLogical} = false       

    % If true, return a list of the solution at each iteration.
    options.return_all {mustBeNumericOrLogical} = false    

    % Initial direction set to search.
    options.direc (:,:) double = eye([length(x0),length(x0)])
end

% Abbreviate variable names.
xtol = options.xtol;
ftol = options.ftol;
maxiter = options.maxiter;
maxfun = options.maxfun;
direc = options.direc;
retall = options.return_all;

% Initialize outputs.
xopt = x0;
fopt = func(x0);
funcalls = 1;
iter = 0;
warnflag = 0;

if retall
    allvecs = {x0};
else
    allvecs = {};
end

% Insure maxiter and maxfun are reasonable.
if maxiter == inf && maxfun == inf
    disp("Warning: maxiter and maxfun are both..." + ...
        " infinite which may cause a loop");
    disp("\n set maxiter to default dimension*1000");
    maxiter = length(x0)*1000;
end

% Outer iteration begins

% Use p0 to record initial point in every iteration.
p0 = x0;  
% Use fopt to record the initial function value in every iteration.
fp0 = fopt; 
% Use init_step to define the init step of linesearch in every iterattion.
deltaf = 0.4*100*max(x0);

while true
% Check the numbers of iterations and function evaluations. 
if funcalls >= maxfun
    warnflag = 1;
    break;
elseif iter >= maxiter 
    warnflag = 2;
    break;
elseif any(isnan(xopt)) || isnan(fopt)
    warnflag = 3;
    break;
end
%%%%%% main optimization process%%%%%%%%%%%%%%%%%%%%%%%%%%%
bigind = 0;
delta = 0;

% Step 1 : Iterate over directions
for i = 1:length(direc)
    d = direc(:, i);
    fx_prev = fopt;
    
    % myfunc = @(x) func(xopt + d * x);
    % [alpha, fopt, ~, feval_count] = direct_search(myfunc, 0, line_step, 0.5, ftol, inf);
    % xopt = xopt + d * alpha;

    [xopt, fopt, feval_count, ~] = linesearch_powell(func, xopt, d, ...
        remaining_calls= (maxfun-funcalls), ...
        init_step=deltaf);

    funcalls = funcalls + feval_count;

    if funcalls >= maxfun
        warnflag = 1;
        return
    end

    % Step 2 : record maximal decents in Step 1.
    if (fx_prev - fopt) > delta
        delta = fx_prev - fopt;
        bigind = i;
    end
end

% record iter times and calculate xopt.
iter = iter + 1;
if retall
    allvecs{end+1} = xopt;
end

% Step3 : Construct extrapolated point, now xopt = pn;
d = xopt - p0; 
x3 = xopt + d;
% Note that x1 = 2pn -p0;
f3 = func(x3); funcalls = funcalls + 1;

if funcalls >= maxfun
    warnflag = 1;
    return
end

% Step 4 : Change directions set or not.
if fp0 > f3
    t = 2 * (fp0 + f3 - 2*fopt);
    temp = (fp0 - fopt - delta);
    t = t * temp^2;
    temp = fp0 - fopt;
    t = t - delta * temp^2;
    if t < 0

        myfunc = @(x) func(xopt + d * x);
        [alpha, fopt, ~, feval_count] = direct_search(myfunc, 0, deltaf, 0.5, ftol, inf);
        xopt = xopt + d * alpha;

        % [xopt, fopt, feval_count, magnitude] = linesearch_powell(...
        %    func, xopt, d, ...
        %    xtol=xtol, ftol = ftol, ...
        %    fval = fopt, remaining_calls= maxfun - funcalls, ...
        %    init_step = deltaf);

        funcalls = funcalls + feval_count;
        if magnitude > 0
            direc(:,bigind) = direc(:,end);
            direc(:, end) = d/norm(d);
        end
    end
end

if funcalls >= maxfun
    warnflag = 1;
    break
end

% Step 5 : Check for convergence
rel_diff = abs(fp0 - fopt) / max([1, abs(fp0), abs(fopt)]);
if rel_diff < ftol
    break;
end

% Steo 6 : Didn't break, set for the next iteration.
deltaf = 0.4*sqrt(abs(fp0-fopt));
p0 = xopt;
fp0 = fopt;

% if all(abs(x1 - xopt) <= 0.1 * xtol) && isnan(terminal_point1)
%     terminal_point1 = xopt;
%     xopt = xopt + 10*xtol;
% else 
%     if all(abs(x1 - xopt) <= 0.1 * xtol) && isnan(terminal_point2)
%         terminal_point2 = xopt;
%         xopt = linesearch_powell(...
%            func, xopt, terminal_point1-terminal_point2 , ...
%            xtol=xtol, ftol = ftol, ...
%            lower_bound = lb, upper_bound = ub, ...
%            fval = fopt, remaining_calls= maxfun - funcalls);
%         if all(abs(terminal_point1 - xopt)<= 0.1 * bnd)...
%                 && all(abs(terminal_point2 - xopt)<= 0.1 * bnd)
%             return
%         else
%             direc(1,:) = terminal_point1 - xopt;
%             terminal_point1 = nan;
%             terminal_point2 = nan;
%         end
%     end
% end


% Final display
if options.disp
    msg = ["Optimization terminated successfully", ...
           "Maximum number of function evaluations exceeded", ...
           "Maximum number of iterations exceeded", ...
           "NaN encountered", ...
           "Result is out of bounds"];
    if warnflag == 0
        status = msg(1);
    else
        status = msg(warnflag + 1);
    end
    fprintf("%s\n", status);
    fprintf("  Final function value: %f\n", fopt);
    fprintf("  Iterations: %d\n", iter);
    fprintf("  Function evaluations: %d\n", funcalls);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

