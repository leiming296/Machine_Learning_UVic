% To apply logistic regression to classify double nonseparable semi-circle
% data. BFGS algorithm is chosen to minimize the logistic loss function.
% Written by W.-S. Lu, University of Victoria. Last modified: Jan. 25, 2015.
% Example:
% [x,y,xp,xn] = data_semi_circle(10,5,-1,1000,9,7);
% [wt,fs] = logistic_BFGS_K(x,y,xp,xn,'logistic_f','logistic_g',zeros(3,1),17,20);
function [wt,fs] = logistic_BFGS_K(x,y,xp,xn,fname,gname,w0,st,K)
N = length(y);
rand('state',st)
N1 = randperm(N);
x0 = x;
y0 = y;
x = x0(:,N1);
y = y0(N1);
d1 = length(w0);
p = [N; d1; y; x(:)];
I = eye(d1);
k = 1;
wk = w0;
Sk = I;
gk = feval(gname,wk,p);
dk = -Sk*gk;
ak = bt_lsearch(wk,dk,fname,gname,p);
dtk = ak*dk;
wk_new = wk + dtk;
while K > k,
      gk_new = feval(gname,wk_new,p);
      gmk = gk_new - gk;
      D = dtk'*gmk;
      if D <= 0,
         Sk = I;
      else
         sg = Sk*gmk;
         sw0 = (1+(gmk'*sg)/D)/D;
         sw1 = dtk*dtk';
         sw2 = sg*dtk';
         Sk = Sk + sw0*sw1 - (sw2'+sw2)/D;
      end
      gk = gk_new;
      wk = wk_new;
      dk = -Sk*gk;
      ak = bt_lsearch(wk,dk,fname,gname,p);
      dtk = ak*dk;
      wk_new = wk + dtk;
      k = k + 1;
end
disp('Final weights:')
wt = wk_new
fs = feval(fname,wt,p);
disp(sprintf('Objective function at solution point: %d.', fs));
Dt = [ones(N,1) x'];
y = y(:);
dwt = (Dt*wt >= 0);
z = dwt + dwt - y - 1;
L = sum(abs(z))/2;
disp(sprintf('Out of N = %d sample points,',N)); 
disp(sprintf('%d points were classified correctly.',N-L));
figure(1)
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
title('Decision Boundary by Logistic Regression (BFGS-BT)')
hold off 