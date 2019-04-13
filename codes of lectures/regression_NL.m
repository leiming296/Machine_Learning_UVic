% To apply linear regression to classify nonseparable double demi-circle
% data using a dilated input space where vetors are dilated to include
% up to 3rd-order features.
% Written by W.-S. Lu, University of Victoria. Last modified: Jan. 25, 2015.
% Example:
% [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% [wt,Ein1,Ein2] = regression_NL(x,y,xp,xn,11);
function [wt,Ein1,Ein2] = regression_NL(x,y,xp,xn,st1)
N = length(y);
rand('state',st1)
N1 = randperm(N);
x0 = x;
y0 = y;
x = x0(:,N1);
y = y0(N1);
z1 = zeros(7,N);
for i = 1:N,
z1(:,i) = [x(1,i)^2; x(1,i)*x(2,i); x(2,i)^2; x(1,i)^3; x(1,i)^2*x(2,i); x(1,i)*x(2,i)^2; x(2,i)^3];
end
z = [x; z1];
Z = [ones(N,1) z'];
y = y(:);
wt = (inv(Z'*Z))*(Z'*y);
disp('Final weights:')
wt
Ein1 = norm(Z*wt-y)^2/N;
disp(sprintf('In-sample error was found to be %d.', Ein1));
dwt = (Z*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
Ein2 = L/N;
disp(sprintf('Out of N = %d sample points,',N));
disp(sprintf('%d points were classified correctly.',N-L));
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