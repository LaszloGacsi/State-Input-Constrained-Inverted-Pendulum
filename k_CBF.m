function outCBF = k_CBF(x)
global par
% outCBF = zeros(2,1);
[h,gradh] = CBF(x);
H = eye(1);
F = [];
Aeq = [];
beq = [];
x0 = 0;
A = [];
B = [];
A = [A; -gradh*g(x); a(x)];
B = [B; par.alpha*h+gradh*f(x)+gradh*g(x)*k_d(x); b(x)];
options = optimoptions('quadprog','Algorithm','active-set');
outCBF = quadprog(H,F,A,B,Aeq,beq,[],[],x0,options) + k_d(x);
end