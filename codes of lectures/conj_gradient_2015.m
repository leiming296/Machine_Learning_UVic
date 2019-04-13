% Conjugate-gradient algorithm for Ax = b with A positive definite.
% Theory: See Practical Optimization Sec. 6.9.
% Input:
% {A,b} -- data for linear equation A and b.
% epsi: termination tolerance.
% Output:
% xs: solution point.
% Written by W.-S. Lu, University of Victoria. Last modified: March 28, 2015.
% ==========================================================
function xs = conj_gradient_2015(A,b,epsi)
n = length(b);
x0 = zeros(n,1);
xk = x0;
rk = b - A*xk;
dk = rk;
rk2_old = rk'*rk;
err = norm(rk);
while  err >= epsi,
    adk = A*dk;
    ak = rk2_old/(dk'*adk);
    xk = xk + ak*dk;
    rk = rk - ak*adk;
    rk2 = rk'*rk;
    bk = rk2/rk2_old;
    rk2_old = rk2;
    dk = rk + bk*dk;
    err = norm(rk);
end
xs = xk;