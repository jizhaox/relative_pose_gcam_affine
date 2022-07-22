%% The affine transformation is calculated from the first-order approximation of the homography
function A = get_affinetransformation(H, x1_p, x2_p)

PointsNum = size(x1_p,2);
A = zeros(2,2,PointsNum);
for i = 1:PointsNum
    
    x1 = x1_p(1,i);
    y1 = x1_p(2,i);
    x2 = x2_p(1,i);
    y2 = x2_p(2,i);
    
    s = H(3,:) * [x1;y1;1];
    
    a11 = (H(1,1) - H(3,1)*x2) / s;
    a12 = (H(1,2) - H(3,2)*x2) / s;
    a21 = (H(2,1) - H(3,1)*y2) / s;
    a22 = (H(2,2) - H(3,2)*y2) / s;
    
    A(:,:,i) = [a11,a12;a21,a22];
end

end

