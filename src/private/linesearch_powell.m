function [xopt, fopt, funcalls, magnitude] = linesearch_powell(func, p, xi, options)
%{
Powell's method to find the minimum of the function ''func(x0 + alpha*direc)''.
Outputs:
xopt : array        -the minimization point.
fopt : number       -the optimal value.
funcalls : int      -numbers of this program evaluate func. 
magnitude : number  -the magnitude of step from init point to xopt.
%}
arguments
    func function_handle        % object function
    p (:,1) double              % initial point
    xi (:,1) double             % search direction

    % Set final accuracy or relative error.
    options.xtol (1,1) double = 1e-6
    options.ftol (1,1) double = 1e-6

    % Init step for search.
    options.init_step (1,1) double = nan
    % Record func(p) to economize one evaluation. 
    options.fval (1,1) double = nan 

    % Get remaining calls for the iteration, it must be at least 1.
    options.remaining_calls (1,1) double {mustBeInteger,...
        mustBePositive} = length(p)* 10000;
end

% Abbreviate variable names.
xi = xi/norm(xi);
myfunc = @(alpha) func(p + alpha * xi);
xtol = options.xtol;      
istep = options.init_step;

f0 = options.fval;
maxfun = options.remaining_calls;
% Initialize outputs.
xopt = p;
if ~isnan(f0)
    fopt = f0;
    funcalls = 0;
else
    fopt = myfunc(0);
    funcalls = 1;
    if funcalls >= maxfun
        return
    end
end
magnitude = 1; % init

% Step 1 : Set default step size bounds
m = 100 * max([abs(xi);1e-12]);

% Step 2 : Gather points for quadratic interpolation.
a = 0;
if isnan(istep)
    b = 0.4*m;
else
    b = istep;
end

fa = fopt;
fb = myfunc(b);funcalls = funcalls + 1;
if funcalls >= maxfun
    return;
end

% Get the third points for interpolation.
if fa > fb
    c = 2*b;
    fc = myfunc(c); funcalls = funcalls + 1;
else
    c = b;fc = fb;
    b = a;fb = fa;
    a = -c;fa = myfunc(a);
    funcalls = funcalls + 1;
end
% Step 3 : Start linesearch process
while funcalls <= maxfun
    n = (a-b)*(c-a)*(b-c);
    nd = fa*(b^2-c^2) + fb*(c^2-a^2) + fc*(a^2-b^2);
    den= fa*(b-c) + fb*(c-a) + fc*(a-b);

    % Check the second derivative n/den, 
    % First If n/den < 0, we finds the minimum
    % or we find the direction to decreasing func to the boundary.
    if (n < 0 && den > 0) || (n > 0 && den <0)
        d = 0.5 * nd/den;
        % Second If d is inside the bounds, we found the minimum
        if abs(d-b) <= m
            % Now we check the convergence criterion
            distances = abs([a, b, c] - d);
            [min_dist, idx] = min(distances);
            % Third If the convergence criterions are met, we will return.
            if min_dist < xtol*1e-3
                if idx == 1
                    xopt = p + a*xi;
                    fopt = fa;
                    magnitude = abs(a);
                elseif idx == 2
                    xopt = p + b*xi;
                    fopt = fb;
                    magnitude = abs(b);
                else
                    xopt = p + c*xi;
                    fopt = fc;
                    magnitude = abs(c);
                end
                return
            % The convergence of x is not met, we continue to
            % iterate, but delete the farthest point a and add d.
            else
                fd = myfunc(d); funcalls = funcalls + 1;
                [a,b,c,fa,fb,fc] = definite_bracket(a,b,c,d,fa,fb,fc,fd);
            end

        % Second Else if d is out of the bounds, we set d to be the 
        % farthest point to decrease func.
        else
            if funcalls >= maxfun
                return;
            end

            if d-b > 0
                d = b+m;
            else
                d = b-m;
            end
                % Use generate function to find the farthest point to d
                % and substitude it by d.
                fd = myfunc(d); funcalls = funcalls + 1;
                [a,b,c,fa,fb,fc] = generate(a,b,c,d,fa,fb,fc,fd);
        end
    % First Else the second derivative is positive, d is a maximal point,
    % we will set d to be the endpoint with the smaller function value.
    else
        if funcalls >= maxfun
            return;
        end
        [f_min, f_max] = quadratic_endpoint_values(a, b, c, fa, ...
            fb, fc, b - m, b + m);
        if f_min <= f_max
            d = b - m;
        else
            d = b + m;
        end
        % Now fd will be used in next line search
        fd = myfunc(d);funcalls = funcalls +1;
        [a,b,c,fa,fb,fc] = generate(a,b,c,d,fa,fb,fc,fd);
    end
end
end





