function [A,gradA] = a(x)
global par
A = zeros(par.N,1);
A(1) = 0;
A(2) = 0;
A(3) = 1;
A(4) = -1;
A(5) = x(2);
A(6) = -x(2);

gradA(1,:) = [0, 0];
gradA(2,:) = [0, 0];
gradA(3,:) = [0, 0];
gradA(4,:) = [0, 0];
gradA(5,:) = [0, 1];
gradA(6,:) = [0, -1];
end