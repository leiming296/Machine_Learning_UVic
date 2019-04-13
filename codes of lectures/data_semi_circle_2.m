% Example: [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% Written by W.-S. Lu, University of Victoria.
% Last modified: Dec. 30, 2014.
function [x,y,xp,xn] = data_semi_circle(r,thk,sep,N,st1,st2)
r1 = r;
r2 = r + thk;
ku = 0;
xn = [];
rand('state',st1)
while ku < N,
xw = 12*randn(2,1);
cr = 0;
xm = norm(xw);
if xm < r2 & xm > r1,
cr = cr + 1;
end
if xw(2) > 0,
cr = cr + 1;
end
if cr == 2,
xn = [xn xw];
ku = ku + 1;
end
end
kd = 0;
xp = [];
rand('state',st2)
while kd < N,
xw = 12*randn(2,1);
cr = 0;
xm = norm(xw);
if xm < r2 & xm > r1,
cr = cr + 1;
end
if xw(2) < 0,
cr = cr + 1;
end
if cr == 2,
xp = [xp xw];
kd = kd + 1;
end
end
ra = 0.5*(r1+r2);
xp(1,:) = xp(1,:) + ra;
xp(2,:) = xp(2,:) - sep;
x = [xp xn];
y = [ones(N,1); -ones(N,1)];
figure(1)
clf
plot(xp(1,:),xp(2,:),'bo')
hold on
plot(xn(1,:),xn(2,:),'r+')
grid
xlabel('\itx_1')
ylabel('\itx_2')
axis([-20 30 -25 25])
axis square
title('Training Data: Double Semi-Circle')
hold off