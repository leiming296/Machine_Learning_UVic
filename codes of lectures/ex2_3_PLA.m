% Perceptron learning algorithm (PLA) as applied to 
% training data in Example 2.3.
% Example: load data_ex2_3
% [wt,t] = ex2_3_PLA(x,y,xp,xn,0.1,[0 0 0]',9);
% Written by W.-S. Lu, University of Victoria.
% Last modified: Jan. 31, 2015.
function [wt,t] = ex2_3_PLA(x,y,xp,xn,r,w0,st)
N = length(y);
t = 0;
wt = w0(:);
D = [ones(N,1) x'];
y = y(:);
dwt = (D*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
if L > 0, 
   flag = 0;
else
   flag = 1;
end
rand('state',st)
while flag == 0,
    Nw = randperm(N);
    tu = 0;
    tv = 0;
    while tu < N & tv == 0,
          ti = Nw(tu + 1);
          di = (D(ti,:)*wt >= 0);
          zi = di + di - y(ti) - 1;
          if zi ~= 0,
             tv = 1;
             t_ind = ti; 
          end
          tu = tu + 1;
    end
    if tv == 0,
       flag = 1;
    else
       wt = wt + (r*y(t_ind))*D(t_ind,:)';  % r is a learning rate
       t = t + 1;
    end
end
dwt = (D*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
disp(sprintf('Solution reached in %d iterations.', t));
disp('Final weights:')
wt
Ein = L/N;
disp(sprintf('In-sample error was found to be %d.', Ein));
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly.',N-L));
figure(1)
subplot(121)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'rx','linewidth',1.5)
grid
xlabel('(a)  \itx_1 (years in residence)')
ylabel('\itx_2 (salary in 10K dollars)')
axis([0 12 0 8])
axis square
title('Training Data')
hold off  
subplot(122)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'rx','linewidth',1.5)
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
title('Decision Boundary Produced by PLA')
hold off  