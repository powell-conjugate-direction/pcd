function [x_opt, f_opt, iter, funcalls] = direct_search(f, x0, alpha0, gamma, tol, max_iter)
% Direct search
% Inpute:
%   f - function handle
%   x0 - initial point
%   alpha0 - initial step
%   gamma - shrink factor (0 < gamma < 1)
%   tol - convergence tolerance
%   max_iter - maximal iteration
% Outputs:
%   x_opt - optimal point
%   f_opt - optimal function value
%   iter - iteration times

% Initialization
x = x0(:); %Ensure x is a row vector.
alpha = alpha0;
n = length(x);
f_current = f(x);
funcalls = 1;
iter = 0;

% Search directions
D = [eye(n), -eye(n)];

while alpha > tol && iter < max_iter
    improved = false;
    
    % Search along all directions
    for k = 1:size(D, 2)
        d = D(:, k);
        x_new = x + alpha * d;
        f_new = f(x_new);
        funcalls = funcalls + 1;
        
        % Check the improvement.
        if f_new < f_current
            x = x_new;
            f_current = f_new;
            improved = true;
            break; 
        end
    end
    
    % No improvements, shrink the step.
    if ~improved
        alpha = gamma * alpha;
    end
    
    iter = iter + 1;
    
    % % Print iteration's info
    % fprintf('Iter %d: f(x) = %.6f, alpha = %.6f\n', iter, f_current, alpha);
end

x_opt = x;
f_opt = f_current;

if alpha <= tol*0.01
    fprintf('Convergence\n');
elseif iter >= max_iter
    fprintf('Stopped: exceed maximal iteration times.\n');
end
end