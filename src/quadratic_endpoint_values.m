function [f_left, f_right] = quadratic_endpoint_values(a, b, c, fa, fb, fc, left, right)
% interpolate 2 quadratic to calculate endpoints vavlue
    
    % Assume f(x) = A*x^2 + B*x + C
    A = ((fa-fb)*(a-c) - (fa-fc)*(a-b))/((a-b)*(a-c)*(b-c));
    B = (fb-fa)/(b-a) - A*(a + b);
    C = fa - A*a^2 - B*a;
    
    % Function values in end point.
    f_left = A*left^2 + B*left + C;
    f_right = A*right^2 + B*right + C;
end