function [B,gradB] = b(x)
global par
B = zeros(par.N,1);
B(1) = 1.75-x(1);
B(2) = 1.75+x(1);
B(3) = par.u_max;
B(4) = -par.u_min;
B(5) = par.P_max;
B(6) = -par.P_min;


gradB(1,:) = [-1, 0];
gradB(2,:) = [1, 0];
gradB(3,:) = [0, 0];
gradB(4,:) = [0, 0];
gradB(5,:) = [0, 0];
gradB(6,:) = [0, 0];
end