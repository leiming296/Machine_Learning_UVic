% For Example 5.1 where linear SVM was applied to the same data as that 
% used in Example 2.3 for PLA.
% Written by W.-S. Lu, University of Victoria.
% Last modified: March 23, 2015.
% Output:
% (w,b): optimal separating weight. 
% sv: support vectors.
% marg: optimal margin of decision boundary.
% Example:
% load data_ex2_3
% [w,b,sv,marg] = svm_linear(x,y,xp,xn);
function [w,b,sv,marg] = svm_linear(x,y,xp,xn)
N = length(y);
ind = 1:1:N;
y = y(:);
d = diag(y)*x';
e = ones(N,1);
cvx_begin quiet
   variable u(N,1);
   dw = d'*u;
   minimize(0.5*(dw'*dw) - e'*u);
   subject to
   y'*u == 0;
   u >= 0;
cvx_end
mu = u;
ind1 = find(mu < 1e-6);
ind2 = setdiff(ind,ind1);
sv = x(:,ind2);
mu1 = mu(ind2);
y1 = y(ind2);
c = mu1(:).*y1;
w = sv*c;
marg = 1/norm(w);
b = y1(1) - w'*sv(:,1);
wt = [b; w];
Dt = [ones(N,1) x'];
dwt = (Dt*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
disp('optimal weights:')
wt
Ein = L/N;
disp(sprintf('In-sample error was found to be %d.', Ein));
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly.',N-L));
disp(sprintf('optimal margin of decision boundary: %d.',marg));
figure(1)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'rx','linewidth',1.5)
nt = size(sv,2);
for i = 1:nt,
    if y1(i) == 1,
        plot(sv(1,i),sv(2,i),'b.','linewidth',1.5)
        plot(sv(1,i),sv(2,i),'b+','linewidth',1.5)
        plot(sv(1,i),sv(2,i),'bx','linewidth',1.5)
    else
        plot(sv(1,i),sv(2,i),'r.','linewidth',1.5)
        plot(sv(1,i),sv(2,i),'r+','linewidth',1.5)
        plot(sv(1,i),sv(2,i),'ro','linewidth',1.5)
    end
end
grid
xlabel('(b)  \itx_1 (years in residence)')
ylabel('\itx_2 (salary in 10K dollars)')
p1 = 0;
p2 = (-wt(2)*p1-wt(1))/wt(3);
q1 = 12;
q2 = (-wt(2)*q1-wt(1))/wt(3);
plot([p1 q1],[p2 q2],'k-','linewidth',1.5)
axis([0 12 0 8])
axis square
title('Decision Boundary Produced by SVM')
hold off  