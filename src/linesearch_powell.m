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
    % Lower bounds for parameter.
    options.lower_bound (:,1) double = -inf(length(p), 1) 
    % upper bounds for parameter.
    options.upper_bound (:,1) double = inf(length(p), 1) 
    % Record func(p) to economize one evaluation. 
    options.fval (1,1) double = nan 
    % Get remaining calls for the iteration, it must be at least 1.
    options.remaining_calls (1,1) double {mustBeInteger,...
        mustBePositive} = length(p)* 1000;
end

% Abbreviate variable names.
xi = xi/norm(xi);
myfunc = @(alpha) func(p + alpha * xi);
xtol = options.xtol;     
ftol = options.ftol;   
istep = options.init_step;
lower = options.lower_bound;
upper = options.upper_bound;
f0 = options.fval;
maxfun = options.remaining_calls;
% Validation the bounds.


% Initialize the points for quadratic interpolation.
% Step 1 : Set problem's bounds and step size bounds

[l_min, l_max] = line_for_search(p, xi, lower, upper);
if isnan(l_min) || isnan(l_max)
    disp("Your bounds are invalid");
    return
end

% Given E, a scalar so that each component of m * |direction_i| <= E * tol
% to define the largest step m in every linesearch.
E = 1e12;
m_com = E * xtol ./ max(abs(xi), 1e-12);  
m = min(m_com);

% Initialize outputs.
xopt = p;
if ~isnan(f0)
    fopt = f0;
    funcalls = 0;
else
    fopt = myfunc(0);
    funcalls = 1;
    if funcalls >= maxfun
        return;
    end
end
magnitude = 1; % init

% Step 2 : Gather points for quadratic interpolation.
a = 0;
if isnan(istep)
    b = min([0.1 * E * xtol, m, l_max]);
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
    c = min([2*b, m, l_max]);
    fc = myfunc(c); funcalls = funcalls + 1;
else
    % Insure that c is inside the maximal step and the given bounds.
    c = b;fc = fb;
    b = a;fb = fa;
    a = max([-b, -m, l_min]); 
    fa = myfunc(a); funcalls = funcalls + 1;
end

% Step 3 : Start linesearch process
max_trials = 200;
for trial = 1:max_trials
    den = (a-b)*(c-a)*(b-c);
    nd = fa*(b^2-c^2) + fb*(c^2-a^2) + fc*(a^2-b^2);
    n= fa*(b-c) + fb*(c-a) + fc*(a-b);

    % Check the second derivative n/den, 
    % First If n/den < 0, we finds the minimum
    % or we find the direction to decreasing func to the boundary.
    if (n < 0 && den > 0) || (n > 0 && den <0)
        d = 1/2 * nd/n;
        % Second If d is inside the bounds, we found the minimum
        if (abs(d-b) <= m) && (l_min <= d) && (d <= l_max)

            % Check the convergence criterion
            distances = abs([a, b, c] - d);
            [min_dist, idx] = min(distances);

            % Third If the convergence criterions are met, we will return.
            if min_dist < xtol*0.05
                if idx == 1
                    xopt = p + a*xi;
                    fopt = fa;
                    magnitude = a;
                elseif idx == 2
                    xopt = p + b*xi;
                    fopt = fb;
                    magnitude = b;
                else
                    xopt = p + c*xi;
                    fopt = fc;
                    magnitude = c;
                end
                return
            % Third Else the convergence criterion of x is not met, we try
            % to verify the convergence criterion of f.
            else
                if funcalls >= maxfun
                    return;
                end
                fd = myfunc(d);funcalls = funcalls + 1;
                % Check remaining calls for func.
                fdistances = abs([fa, fb, fc] - fd);
                [min_dist, idx] = min(fdistances);
                if min_dist < ftol * 0.03
                    if idx == 1
                        xopt = p + a*xi;
                        fopt = fa;
                        magnitude = a;
                    elseif idx == 2
                        xopt = p + b*xi;
                        fopt = fb;
                        magnitude = b;
                    else
                        xopt = p + c*xi;
                        fopt = fc;
                        magnitude = c;
                    end
                    return
                end
                % The convergence of x and f are not met, we continue to
                % iterate.
                [a,b,c,fa,fb,fc] = generate(a,b,c,d,fa,fb,fc,fd);
            end

        % Second Else if d is out of the bounds, we set d to be the 
        % farthest point to decrease func.
        else
            if funcalls >= maxfun
                return;
            end

            if d-a > 0 
                d = min([l_max, a+m]);
            else
                d = max([a-m, l_min]);
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
        interval_min = max(a - m, l_min);
        interval_max = min(a + m, l_max);
        [f_min, f_max] = quadratic_endpoint_values(a, b, c, fa, ...
            fb, fc, interval_min, interval_max);
        if f_min <= f_max
            d = interval_min;
        else
            d = interval_max;
        end
        % Now fd will be used in next line search
        fd = myfunc(d);funcalls = funcalls +1;
        [a,b,c,fa,fb,fc] = generate(a,b,c,d,fa,fb,fc,fd);
    end
end
end





