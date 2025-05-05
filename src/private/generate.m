function [na,nb,nc,fna,fnb,fnc] = generate(a,b,c,d,fa,fb,fc,fd)
% To find the farthest point to d and substitude it by d.
na = a;fna = fa;
nb = b;fnb = fb;
nc = c;fnc = fc;
% We need to delete the point farthest to d and add d. 
distances = abs([a, b, c] - d);
[~, idx] = max(distances);
if idx == 1
    na = d;
    fna = fd;
elseif idx == 2
    nb = d;
    fnb = fd;
else
    nc = d;
    fnc = fd;
end
[x,idx] = sort([na,nb,nc]);
f = [fna,fnb,fnc];
fs = f(idx);
na = x(1);nb =x(2);nc=x(3);
fna= fs(1);fnb = fs(2);fnc = fs(3);
end