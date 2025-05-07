format long g
% f = @(x) exp(-x) + x^2;
% f = @(x) sin(x) + 0.1*x^2;
% f = @(x) abs(x - 3);
f = @(x) (1 - x).^2 + 10*(x - x.^2).^2;
% Set parameters;
x0 = 4;


alpha0 = 0.1;      % initial step
gamma = 0.5;       % shrink factor
tol = 1e-6;        % tolerance
max_iter = 5000;   % maximal iteration times

[x_opt1, f_opt1,funcalls1,~,allb,alla,allc,alld] = linesearch_powell(f, x0, 1);


x = linspace(-10,10,1000);
y = f(x);
hold on
plot(x,y,'-');

yb = f(x0 + allb);ya = f(alla+x0);yc = f(allc+x0);
plot(x_opt1 + x0,f_opt1, '*',MarkerSize = 10);

% semilogy(alla,ya,'o',Color='r',MarkerSize=10)
plot(allb,yb,'.',Color='r',MarkerSize=10);
% semilogy(allc,yc,'^',Color='r',MarkerSize=10);

color = ['r','b','g','cyan','black'];

% 输入三个点
for i =1:5
    points = [alla(i), ya(i); allb(i), yb(i); allc(i), yc(i)];
    a = points(1,1);b = points(2,1);c =points(3,1);
    % disp([a,b,c]);
    fa = points(1,2);fb = points(2,2);fc = points(3,2);
    % nd = fa*(b^2-c^2) + fb*(c^2-a^2) + fc*(a^2-b^2);
    % den= fa*(b-c) + fb*(c-a) + fc*(a-b);

    % 计算二次插值系数 ax² + bx + c
    A = [points(:,1).^2, points(:,1), ones(3,1)]; % 构建系数矩阵
    coeffs = A \ points(:,2); % 解线性方程组
    a = coeffs(1); b = coeffs(2); c = coeffs(3);
    g = @(x) a*x.^2+b*x + c;
    plot(x,g(x),'-',Color=color(i));

end
% [x_opt, f_opt,iter, funcalls] = direct_search(f, x0, alpha0, gamma, tol, max_iter)
fprintf("optimal point :%6f\n", x_opt1);
fprintf("optimal value :%6f", f_opt1);
fprintf("function calls: %6f", funcalls1);
