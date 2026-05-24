function [hb,gradhb] = h_b(x)
global par
hb = par.c - x(1)^2 - (1/(2*par.mu))*(x(2)+par.K_b*x(1))^2;
gradhb = [-2*x(1) - (1/par.mu)*(x(2)+par.K_b*x(1))*par.K_b, -(1/par.mu)*(x(2)+par.K_b*x(1))];
end