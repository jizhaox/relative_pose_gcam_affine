function [PT, Distance, Nvector] = generate_3Dscenepoints(n_point)
% Input
%  n_point: the number of ACs, i.e., the number of 3D plane
% Output
%  PT: (1+4)*n_point 3D scene points. At each 3D plane, we generate 1 PT (PC of AC) + 4 PTs for computing the affine transformation
%  D: the distance between the 3D plane and the body reference at time i
%  Nvector: the normal vector at the body reference at time i

PT = zeros(3,5,n_point);
Distance = zeros(1,n_point);
Nvector = zeros(3,n_point);

for i = 1:n_point
    Distance(i) = rand(1);
    Nvector(:,i) = randn(3,1);
    Nvector(:,i) = Nvector(:,i)/norm(Nvector(:,i));
    % set Y and Z coordinates of PTs 
    PT(2:3,:,i) = rand(2,5);
    
    for j = 1:size(PT,2)
        PT(1,j,i) = (Distance(i) - PT(2,j,i)*Nvector(2,i) - PT(3,j,i)*Nvector(3,i))/Nvector(1,i);
    end
end
