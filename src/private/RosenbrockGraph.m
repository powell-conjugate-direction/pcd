% 定义计算范围
x = linspace(-2, 2, 100);
y = linspace(-1, 3, 100);
[X, Y] = meshgrid(x, y);

% 计算函数值
Z = (1 - X).^2 + 100*(Y - X.^2).^2;

% 绘制3D曲面图
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Rosenbrock Function (3D Surface)');
colormap('jet'); colorbar;
view(-20, 30); % 调整视角

% 绘制等高线图
figure;
contour(X, Y, Z, logspace(-1, 3, 50));
xlabel('x'); ylabel('y');
title('Rosenbrock Function (Contour)');
hold on;
plot(1, 1, 'rp', 'MarkerSize', 10); % 标记最小值点
plot(-1,1,'*','MarkerSize',10);
quiver(-1,1, 1,0, 0.5, 'k', 'LineWidth',1.5, 'MaxHeadSize',0.5); % 指向(-1,1)
text(1.1, 1.1, 'Global Minimum (1,1)', 'Color', 'red');
grid on;