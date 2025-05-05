% Test Rosenbrock
rosenbrock = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

% Set parameters
x0 = [-1.2; 1];    % inital point
alpha0 = 0.1;      % initial step
gamma = 0.5;       % shrink factor
tol = 1e-6;        % tolerance
max_iter = 5000;   % maximal iteration times

[x_opt, f_opt, iter] = direct_search(rosenbrock, x0, alpha0, gamma, tol, max_iter);

% Print results.
fprintf('最优解: x = [%.6f, %.6f]\n', x_opt(1), x_opt(2));
fprintf('最优函数值: %.6f\n', f_opt);
fprintf('迭代次数: %d\n', iter);