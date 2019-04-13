%  Program: inex_lsearch.m
%  Title: Inexact Line Search
%  Description: Implements Fletcher's inexact line search described in 
%  Algorithm 4.6.
%  Theory: See Practical Optimization Sec. 4.8
%  Input:
%    x:  initial point
%    s:  search direction
%    F:  objective function to be minimized along the direction of s  
%    G:  gradient of objective function F 
%   p1:  internal parameters that are required for the implementation of
%        the line search regardless of the application at hand.
%        It is a string (e.g. 'rho = 0.1') and can be a combination several
%        internal parameters (e.g., 'rho = 0.25; sigma=0.5').
%        Useful p1's include:                            default value
%        'rho=   ' defines right bracket                      0.1
%        'sigma= ' defines left bracket  (sigma >= rho)       0.1
%        'tau= '   defines minimum step for sectioning        0.1
%        'chi='                                               0.75
%   p2:  user-defined parameter vector. Note that p2 must be a vector
%        with all components numerically specified. The order in which
%        the components of p2 appear must be the same as what they appear 
%        in function F and gradient G. For example, if p2 = [a b], then 
%        F.m and G.m must be in the form of function z = F(x,p2)and
%        function z = G(x,p2)
%  Output:
%     z: acceptable value of alpha
%  Example 1: 
%  Perform inexact line search using the Himmelblau function
%     f(x1,x2) = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%  starting from point xk = [6 6]' along search direction s = [-1 -1]'
%  using the default parameter values. 
%  Solution: 
%  Execute the command
%     a1 = inex_lsearch([6 6]',[-1 -1]','f_himm','g_himm')
%  Example 2: 
%  Perform inexact line search using the Himmelblau funtion
%  starting from point xk = [6 6]' along search direction s = [-1 -1]'
%  with sigma = 0.5.
%  Solution: 
%  Execute the command
%     a2 = inex_lsearch([6 6]',[-1 -1]','f_himm','g_himm','sigma = 0.5')
%  Example 3: 
%  Perform inexact line search using the paramerized Himmelblau funtion
%     f(x1,x2,a,b) = (x1^2 + x2 - a^2)^2 + (x1 + x2^2 - b^2)^2
%  starting from point xk = [6 6]' along search direction s = [-1 -1]',
%  with parameters a = 3.2 and b = 2.6.
%  Solution: 
%  Execute the command
%     a3 = inex_lsearch([6 6]',[-1 -1]','f_himm_p','g_himm_p',[3.2 2.6])
%  Notes:
%  1. Command 
%       z = inex_lsearch(xk,s,F,G,p1,p2)
%     adds a new function inex_lsearch to MATLAB's vocabulary.
%  2. Do not use a semicolon in commands 
%      a1 = inex_lsearch([6 6]',[-1 -1]','f_himm','g_himm')
%      a2 = inex_lsearch([6 6]',[-1 -1]','f_himm','g_himm','sigma = 0.5')
%      a3 = inex_lsearch([6 6]',[-1 -1]','f_himm_p','g_himm_p',[3.2 2.6])
%    otherwise the acceptable value of alpha will not be displayed.
% 3. f_himm and g_himm are user defined functions implemented in m-files  
%    f_himm.m and g_himm.m, respectively, and are used to evaluate the 
%    Himmelblau function and its derivative. Similarly, f_himm_p and 
%    g_himm_p are user defined functions for the parameterized Himmelblau
%    function and its gradient.
% 4. If you plan to use this line search as part of an optimization 
%    algorithm delete lines 73, 74, and 171.
%==========================================================================
function z = inex_lsearch(xk,s,F,G,p1,p2)
k = 0;
m = 0;
tau = 0.1;
chi = 0.75;
rho = 0.1;
sigma = 0.1;
mhat = 400;
epsilon = 1e-10;
xk = xk(:);
s = s(:);
parameterstring ='';
% evaluate given parameters:
  if nargin > 4,
    if isstr(p1),
      eval([p1 ';']);
    else
      parameterstring = ',p1';
    end
  end
  if nargin > 5,
    if isstr(p2),
      eval([p2 ';']);
    else
      parameterstring = ',p2';
    end
  end
% compute f0 and g0
  eval(['f0 = ' F '(xk' parameterstring ');']);
  eval(['gk = ' G '(xk' parameterstring ');']);
  m = m+2;
  deltaf0 = f0;
% step 2 Initialize line search
  dk = s;
  aL = 0;
  aU = 1e99;
  fL = f0;
  dfL = gk'*dk;
  if abs(dfL) > epsilon,
    a0 = -2*deltaf0/dfL;
  else
    a0 = 1;
 end
 if ((a0 <= 1e-9)|(a0 > 1)),
    a0 = 1;
 end
%step 3
 while 1,
    deltak = a0*dk;
    eval(['f0 = ' F '(xk+deltak' parameterstring ');']);
    m = m + 1;
%step 4
    if ((f0 > (fL + rho*(a0 - aL)*dfL)) & (abs(fL - f0) > epsilon) & (m < mhat))
      if (a0 < aU)
        aU = a0;
      end
      % compute a0hat using equation 7.65
      a0hat = aL + ((a0 - aL)^2*dfL)/(2*(fL - f0 + (a0 - aL)*dfL));
      a0Lhat = aL + tau*(aU - aL);
      if (a0hat < a0Lhat)
        a0hat = a0Lhat;
      end
      a0Uhat = aU - tau*(aU - aL);
      if (a0hat > a0Uhat)
        a0hat = a0Uhat;
      end
      a0 = a0hat;
    else
      eval(['gtemp =' G '(xk+a0*dk' parameterstring ');']);
      df0 = gtemp'*dk;
      m = m + 1;
      % step 6
      if (((df0 < sigma*dfL) & (abs(fL - f0) > epsilon) & (m < mhat) & (dfL ~= df0)))
        deltaa0 = (a0 - aL)*df0/(dfL - df0);
        if (deltaa0 <= 0)
          a0hat = 2*a0;
        else
          a0hat = a0 + deltaa0;
        end
        a0Uhat = a0 + chi*(aU - a0);
        if (a0hat > a0Uhat)
          a0hat = a0Uhat;
        end
        aL = a0;
        a0 = a0hat;
        fL = f0;
        dfL = df0;
      else
        break;
      end
    end
 end % while 1
 if a0 < 1e-5,
    z = 1e-5;
 else
    z = a0;
 end 