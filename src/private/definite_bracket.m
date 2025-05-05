function [na,nb,nc,fna,fnb,fnc,cand] = definite_bracket(a,b,c,d,fa,fb,fc,fd)
% We delete one point among {a,b,c} to get a definite bracket. If there is
% no definite bracket, we will delete the point with the largest function
% value.
x  = [a, b, c, d];
f = [fa, fb, fc, fd];
[~, idx] = sort(x);
fs = f(idx);
cand = [];

if (fs(2) > fs(3)) && (fs(3) < fs(4))
    cand = [cand, idx(1)];
end

if (fs(1) > fs(3)) && (fs(3) < fs(4))
    cand = [cand, idx(2)];
end

if (fs(1) > fs(2)) && (fs(2) < fs(4))
    cand = [cand, idx(3)];
end

if (fs(1) > fs(2)) && (fs(2) < fs(3))
    cand = [cand, idx(4)];
end

cand = setdiff(cand,4);

if ~isempty(cand)
    [~, mi] = max(f(cand));
    di = cand(mi);
else
    [~, di] = max([fa, fb, fc]);
end

% delete the point
x(di) = [];
f(di) = [];

% outputs
na = x(1);fna = f(1);
nb = x(2);fnb = f(2);
nc = x(3);fnc = f(3);
end