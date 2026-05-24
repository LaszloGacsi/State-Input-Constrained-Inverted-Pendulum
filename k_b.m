function [KB, gradKB] = k_b(x)
global par
KB = -sin(x(1))-par.K*x(2);
gradKB = [-cos(x(1)), -par.K];
end
