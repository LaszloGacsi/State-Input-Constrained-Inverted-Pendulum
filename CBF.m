function [h,gradh] = CBF(x)
global par
[KB, gradKB] = k_b(x);
[A, gradA] = a(x);
[B,gradB] = b(x);
h = zeros(par.N,1);
gradh = zeros(par.N,par.n);
    for i = 1:par.N
        h(i) = B(i) - A(i)*KB;
        gradh(i,:) = gradB(i,:) - gradA(i,:)*KB - A(i)*gradKB;
    end
end