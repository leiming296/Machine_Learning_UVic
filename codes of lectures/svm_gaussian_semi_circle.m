% For Example 5.3 where SVM with a Gaussian kernel was applied to the same 
% data as that produced in Example 3.1.
% Written by W.-S. Lu, University of Victoria. Last modified: March 26, 2015.
% Input:
% (x,y); training data.
% xp: training data with positive labels.
% xn: training data with negative labels.
% st: initial random state for mixing training data.
% sig: \sigma for Gaussian kernel.
% Output:
% mu: optimal solution of the dual problem.
% b: optimal intersect. 
% sv: "support vectors" in basic feature space.
% cpt: CPU time consumed by CVX.
% Example:
% [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% [mu,b,sv] = svm_gaussian_semi_circle(x,y,xp,xn,17,0.8);
function [mu,b,sv,cpt] = svm_gaussian_semi_circle(x,y,xp,xn,st,sig)
N = length(y);
rand('state',st)
N1 = randperm(N);
x0 = x;
y0 = y;
x = x0(:,N1);
y = y0(N1);
y = y(:);
D1 = eye(N);
for i = 1:(N-1),
    for j = (i+1):N,
        nij = (norm(x(:,i)-x(:,j)))^2;
        D1(i,j) = exp(-nij/(2*sig^2));
        D1(j,i) = D1(i,j);
    end
end
Y = y*y';
D = Y.*D1;
clear Y 
e = ones(N,1);
cvx_begin quiet
   variable u(N,1);
   minimize(0.5*u'*D*u - e'*u);
   subject to
   u'*y == 0;
   u >= 0;
cvx_end
cpt = cvx_cputime
mu = u;
ind = 1:1:N;
ind1 = find(mu < 1e-5);
ind2 = setdiff(ind,ind1);
sv = x(:,ind2);
mu1 = mu(ind2);
y1 = y(ind2);
c = mu1(:).*y1;
di = D1(ind2(1),ind2);
b = y1(1) - di*c;
nt = length(ind2);
dw = b*ones(N,1);
for i = 1:N,
    xi = x(:,i);
    dwi = dw(i);
    for j = 1:nt,
        nij = (norm(sv(:,j)-xi))^2;
        dij = exp(-nij/(2*sig^2));
        dwi = dwi + c(j)*dij;
    end
    dw(i) = dwi;
end
dwt = (dw >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
Ein = L/N;
disp(sprintf('In-sample error was found to be %d.', Ein));
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly.',N-L));
N2 = 100;
[x1,x2] = meshgrid(-20:50/N2:30,-25:50/N2:25);
h = b*ones(N2+1,N2+1);
for i = 1:nt,
    x1w = sv(1,i);
    x2w = sv(2,i);
    m1w = (x1w - x1).^2;
    m2w = (x2w - x2).^2;
    m2 = sqrt(m1w + m2w)/(2*sig^2);
    Kw = exp(-m2);
    h = h + c(i)*Kw;
end
figure(1)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'r+','linewidth',1.5)
v = -1e-6:1e-6:1e-6;
contour(x1,x2,h,v,'k-','linewidth',1.5);
grid
xlabel('\itx_1')
ylabel('\itx_2')
axis square
axis([-20 30 -25 25])
hold off