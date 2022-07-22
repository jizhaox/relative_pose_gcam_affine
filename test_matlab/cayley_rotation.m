function R = cayley_rotation(cay)

R = zeros(3, 3);
R(1,1) = 1 + cay(1)^2 - cay(2)^2 - cay(3)^2 ;
R(2,2) = 1 - cay(1)^2 + cay(2)^2 - cay(3)^2 ;
R(3,3) = 1 - cay(1)^2 - cay(2)^2 + cay(3)^2 ;
R(1,2) = 2 * (cay(1)*cay(2) - cay(3)) ;
R(1,3) = 2 * (cay(1)*cay(3) + cay(2)) ;
R(2,1) = 2 * (cay(1)*cay(2) + cay(3)) ;
R(2,3) = 2 * (cay(2)*cay(3) - cay(1)) ;
R(3,1) = 2 * (cay(1)*cay(3) - cay(2)) ;
R(3,2) = 2 * (cay(2)*cay(3) + cay(1)) ;

R = R/(1 + cay'*cay);