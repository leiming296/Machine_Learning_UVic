%  Program: bt_lsearch.m
%  Title: Inexact Line Search by Backtracking
%  Description: Implements inexact line search described in 
%  Sec. 9.2 of Boyd-Vanderberghe's book.
%  Theory: Sec. 9.2 of Boyd-Vanderberghe's book.
%  Input:
%    x:  initial point
%    s:  search direction
%    F:  objective function to be minimized along the direction of s  
%    G:  gradient of the objective function.
%   p1:  internal parameters that are required for the implementation of
%        the line search regardless of the application at hand.
%        It is a string (e.g. 'rho = 0.1') and can be a combination several
%        internal parameters (e.g., 'rho = 0.2; gma = 0.6').
%        Useful p1's include:                            default value
%        'rho=   ' defines                                    0.1
%        'gma=   ' defines                                    0.5
%   p2:  user-defined parameter vector. Note that p2 must be a vector
%        with all components numerically specified. The order in which
%        the components of p2 appear must be the same as what they appear 
%        in function F. For example, if p2 = [a b], then F.m must be in the 
%        form of function z = F(x,p2).
%  Output:
%    a:  acceptable value of alpha
%  Example 1: 
%  Perform inexact line search using the Himmelblau function
%     f(x1,x2) = (x1^2 + x2 - 11)^2 + (x1 + x2^2 - 7)^2
%  starting from point xk = [6 6]' along search direction s = [-1 -1]'
%  using the default parameter values. 
%  Solution: 
%  Execute the command
%     a1 = bt_lsearch([6 6]',[-1 -1]','f_himm','g_himm')
%  Example 2: 
%  Perform inexact line search using the Himmelblau funtion
%  starting from point xk = [6 6]' along search direction s = [-1 -1]'
%  with rho = 0.4.
%  Solution: 
%  Execute the command
%     a2 = bt_lsearch([6 6]',[-1 -1]','f_himm','g_himm','rho = 0.4')
%  Example 3: 
%  Perform inexact line search using the paramerized Himmelblau funtion
%     f(x1,x2,a,b) = (x1^2 + x2 - a^2)^2 + (x1 + x2^2 - b^2)^2
%  starting from point xk = [6 6]' along search direction s = [-1 -1]',
%  with parameters a = 3.2 and b = 2.6.
%  Solution: 
%  Execute the command
%     a3 = bt_lsearch([6 6]',[-1 -1]','f_himm_p','g_himm_p',[3.2 2.6])
%==========================================================================
function a = bt_lsearch(x,s,F,G,p1,p2)
rho = 0.1;
gma = 0.5;
x = x(:);
s = s(:);
a = 1;
parameterstring ='';
% evaluate given parameters:
  if nargin > 4,
    if ischar(p1),
      eval([p1 ';']);
    else
      parameterstring = ',p1';
    end
  end
  if nargin > 5,
    if ischar(p2),
      eval([p2 ';']);
    else
      parameterstring = ',p2';
    end
  end
 eval(['f0 = ' F '(x' parameterstring ');']);
 eval(['g0 = ' G '(x' parameterstring ');']);
 eval(['f1 = ' F '(x+a*s' parameterstring ');']);
 f2 = f0 + rho*a*g0'*s;
 er = f1 - f2;
 while er > 0,
     a = gma*a;
     eval(['f1 = ' F '(x+a*s' parameterstring ');']);
     f2 = f0 + rho*a*g0'*s;
     er = f1 - f2;
 end
 if a < 1e-5,
    a = 1e-5;
 end 