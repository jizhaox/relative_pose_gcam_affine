function [R, t] = find_optimal_transformation(s1, s2)
% Input
%  s1 and s2: focal points of the perspective cameras
%  in the original frame
% Output
%  R and t satisfy: the direction of R*s1+t and R*s2+t
%  are [-1;-1;-1] and [1;1;1], respectively.
t = (s1+s2)/2;
a = s2-s1;
a = a/norm(a);
b = [1;1;1]/sqrt(3);
v = cross(a,b);
s = norm(v);
v = v/s;
if (s<1e-10)
    R = eye(3);
    t = -R*t;
    return;
end
SV = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = eye(3) + SV*s + SV^2*(1-dot(a,b));
t = -R*t;