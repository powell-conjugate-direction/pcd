function [lmin, lmax] = line_for_search(x0, alpha, lower_bound, upper_bound)
%{
To get the range of step l for initial point x0, such that
 x0 + alpha * l in the given bound
outputs:
   lmin: least step 
   lmax: maximal step
%}
arguments
    x0 (:,1) double             % initial point
    alpha (:,1) double          % search direction
    lower_bound (:,1) double    % lower bound for every components
    upper_bound (:,1) double    % upper bound for every components
end

% check whether the direction is zero.
if all(alpha == 0) 
    lmin = 0;
    lmax = 0;
    return
end

% select nonzero components in the direction.
nonzero = (alpha ~= 0);
alpha_nz = alpha(nonzero);
x0_nz = x0(nonzero);
lb_nz = lower_bound(nonzero);
ub_nz = upper_bound(nonzero);

% calculate the magnitude of bounds for every nonzero components 
low = (lb_nz - x0_nz) ./ alpha_nz;
high = (ub_nz - x0_nz) ./ alpha_nz;

% check the lmin and lmax with sign of components of alpha.
pos = alpha_nz > 0;
low_pos = low;
low_neg = low;
high_pos = high;
high_neg = high;
low_pos(pos == 0) = 0;
low_neg(pos == 1) = 0;
high_pos(pos == 0) = 0;
high_neg(pos == 1) = 0;

lmin = low_pos  + high_neg;
lmax = low_neg + high_pos;

% get the real bounds.
lmin = max(lmin);
lmax = min(lmax);

% Validate the result. 
if isfinite(lmax) && isfinite(lmin) && lmax < lmin
    lmin = nan;
    lmax = nan;
    return
elseif lmax == -inf || lmin == inf
    lmin = nan;
    lmax = nan;
    return
end

end