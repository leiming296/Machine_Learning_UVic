% For Example 5.2 where SVM with 3rd-order nonlinearity was
% applied to the same data as that produced in Example 3.1.
% Written by W.-S. Lu, University of Victoria.
% Last modified: March 23, 2015.
% Output:
% (w,b): optimal separating weight. 
% sv: "support vectors" in basic feature space.
% Example:
% [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% [w,b,sv] = svm_NL_semi_circle(x,y,xp,xn,17);
function [w,b,sv] = svm_NL_semi_circle(x,y,xp,xn,st)
N = length(y);
rand('state',st)
N1 = randperm(N);
x0 = x;
y0 = y;
x = x0(:,N1);
y = y0(N1);
y = y(:);
z1 = zeros(7,N);
for i = 1:N,
    z1(:,i) = [x(1,i)^2; x(1,i)*x(2,i); x(2,i)^2; x(1,i)^3; x(1,i)^2*x(2,i); x(1,i)*x(2,i)^2; x(2,i)^3];
end
z = [x; z1];
D = diag(y)*z';
e = ones(N,1);
cvx_begin quiet
   variable w(9,1);
   variable b(1)
   minimize(w'*w);
   subject to
   D*w + y*b >= e;
cvx_end
c = D*w + y*b;
ind1 = find(abs(c - e) < (1e-5)*e);
sv = x(:,ind1);
y1 = y(ind1);
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
[x1,x2] = meshgrid(-20:50/100:30,-25:50/100:25);
w0 = wt(1); w1 = wt(2); w2 = wt(3); w3 = wt(4); w4 = wt(5); w5 = wt(6); 
w6 = wt(7); w7 = wt(8); w8 = wt(9); w9 = wt(10);
h1 = w0 + w1*x1 + w2*x2 + w3*(x1.^2) + w4*(x1.*x2) + w5*(x2.^2);
h = h1 + w6*(x1.^3) + w7*((x1.^2).*x2) + w8*(x1.*(x2.^2)) + w9*(x2.^3);
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