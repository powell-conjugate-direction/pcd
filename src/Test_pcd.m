% Test on Rosenbrock
rosenbrock = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
x0 = [-1.2,1];

[x_opt, f_opt, direc, iter, warnflag, funcalls, allvecs] = ...
    powell_conjugate_direction(rosenbrock, x0, return_all=true)
