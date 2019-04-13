% For Example 5.4 where soft-margin SVM with 2rd-order nonlinearity
% was applied to the same data as that produced in Example 4.6.
% Written by W.-S. Lu, University of Victoria. Last modified: March 27, 2015.
% Input:
% (x,y): training data.
% tau: parameter \tau for the total violation term.
% Output:
% (w,b): optimal separating weight. 
% Example:
% load data_ex5_4
% [w,b] = svm_soft_NL(x,y,xp,xn,0.7);
function [w,b] = svm_soft_NL(x,y,xp,xn,tau)
N = length(y);
y = y(:);
z1 = zeros(3,N);
for i = 1:N,
    z1(:,i) = [x(1,i)^2; x(1,i)*x(2,i); x(2,i)^2];
end
z = [x; z1];
D = diag(y)*z';
e = ones(N,1);
cvx_begin quiet
   variable w(5,1);
   variable b(1);
   variable s(N,1);
   minimize(0.5*w'*w + tau*sum(s));
   subject to
   D*w + y*b + s >= e;
   s >= 0;
cvx_end
wt = [b; w];
Dt = [ones(N,1) z'];
dwt = (Dt*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
disp('optimal weights:')
wt
Ein = L/N;
disp(sprintf('In-sample error was found to be %d.', Ein));
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly.',N-L));
figure(1)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'rx','linewidth',1.5)
[x1,x2] = meshgrid(0:1/100:1,0:1/100:1);
w0 = wt(1); w1 = wt(2); w2 = wt(3); w3 = wt(4); w4 = wt(5); w5 = wt(6);
h = w0 + w1*x1 + w2*x2 + w3*(x1.^2) + w4*(x1.*x2) + w5*(x2.^2);
figure(1)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'r+','linewidth',1.5)
v = -1e-6:1e-6:1e-6;
contour(x1,x2,h,v,'k-','linewidth',1.5);
grid
xlabel('Average intensity')
ylabel('Symmetry')
axis square
axis([0 1 0 1])
hold off