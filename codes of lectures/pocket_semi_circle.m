% Pocket algorithm (PCA) as applied to training data known as semi-circle.
% Example:
% [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% [wt,t,ein_pocket] = pocket_semi_circle(x,y,xp,xn,0.1,[0 0 0]',15,26,40);
% Written by W.-S. Lu, University of Victoria. % Last modified: Jan. 25, 2015.
function [wt,t,ein_pocket] = pocket_semi_circle(x,y,xp,xn,lr,w0,st1,st2,K)
N = length(y);
rand('state',st1)
N1 = randperm(N);
x0 = x;
y0 = y;
x = x0(:,N1);
y = y0(N1);
t = 0;
w = w0(:);
w_pocket = w;
D = [ones(N,1) x'];
y = y(:);
dwt = (D*w >= 0);
acc = sum(y == 2*dwt - 1);
acc_pocket(1) = acc;
rand('state',st2)
N1 = randperm(N);
while t < K & acc_pocket(t+1) < N,
     ii = mod(t,N) + 1;
     i = N1(ii);
     fi = 2*dwt(i) - 1;
     w = w + lr*(y(i)-fi)*D(i,:)';
     dwt = (D*w >= 0);
     acc = sum(y == 2*dwt - 1);
     if acc > acc_pocket(t+1),
        w_pocket = w;
        acc_pocket(t+2) = acc;
     else
        acc_pocket(t+2) = acc_pocket(t+1);
     end
     t = t + 1;
end
disp(sprintf('Solution reached in %d iterations.', t));
disp('Final weights:')
wt = w_pocket
ent = N - acc_pocket(end);
Ein = ent/N;
disp(sprintf('In-sample error was found to be %d.', Ein));
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly',N - ent));
figure(1)
subplot(121)
plot(xp(1,:),xp(2,:),'bo','linewidth',1.5)
hold on
plot(xn(1,:),xn(2,:),'r+','linewidth',1.5)
grid
xlabel('\itx_1')
ylabel('\itx_2')
p1 = -20;
p2 = (-wt(2)*p1-wt(1))/wt(3);
q1 = 30;
q2 = (-wt(2)*q1-wt(1))/wt(3);
plot([p1 q1],[p2 q2],'k-','linewidth',1.5)
axis([-20 30 -25 25])
axis square
title('Decision Boundary Produced by PCA')
hold off
subplot(122)
ein_pocket = (N - acc_pocket)/N;
plot(1:1:t+1,ein_pocket,'b-','linewidth',1.5)
axis square
axis([1 t+1 0 1.1*ein_pocket(1)])
xlabel('number of iterations')
ylabel('in-sample error')
grid