function [x,xhist] = LBFGSB(func,x0,l,u,options)
% function [x,xhist] = LBFGSB(func,x0,l,u,options)
% Perform bound-constrained optimization with L-BFGS-B.
% INPUTS:
%  x0: [n,1] initial design vector.
%  l: [n,1] lower bound constraint vector.
%  u: [n,1] upper bound constraint vector.
%  options: matlab struct with the optional components:
%   'm': the maximum number of stored L-BFGS iteration pairs.
%   'tol': the convergence tolerance for the projected gradient.
%   'display': true/false - should iterations be displayed?
%   'xhist': true/false - should the entire search history be stored?

% validate the inputs
[x0] = validate_inputs(func,x0,l,u);

% set options
[m,tol,max_iters,display,xhistory] = set_options(options);

% initialize BFGS variables
n = length(x0);
Y = zeros(n,0);
S = zeros(n,0);
W = zeros(n,1);
M = zeros(1,1);
theta = 1;

% initialize objective variables
x = x0;
[f,g] = feval(func, x);

% initialize quasi-Newton iterations
k = 0;

% print out some useful information, if specified
if (display)
  fprintf(' iter        f(x)          optimality\n')
  fprintf('-------------------------------------\n')
  opt = get_optimality(x,g,l,u);
  fprintf('%3d %16.8f %16.8f\n',k,f,opt);
end

% save the xhistory, if specified
xhist = [];
if (xhistory)
  xhist = [xhist x0];
end

% perform quasi-Newton iterations
while ( (get_optimality(x,g,l,u) > tol) && (k < max_iters) )
  
  % update search information
  x_old = x;
  g_old = g;
  
  % compute the new search direction
  [xc, c] = get_cauchy_point(x,g,l,u,theta,W,M);
  [xbar, line_search_flag] = subspace_min(x,g,l,u,xc,c,theta,W,M);

  alpha = 1.0;
  if (line_search_flag)
    [alpha] = strong_wolfe(func,x,f,g,xbar-x);
  end
  x = x + alpha * (xbar - x);
  
  % update the LBFGS data structures
  [f,g] = feval(func, x);
  y = g - g_old;
  s = x - x_old;
  curv = abs(transpose(s)*y);
  if (curv < eps)
    fprintf(' warning: negative curvature detected\n');
    fprintf('          skipping L-BFGS update\n');
    k = k+1;
    continue;
  end
  if (k < m)
    Y = [Y y];
    S = [S s];
  else
    Y(:,1:m-1) = Y(:,2:end);
    S(:,1:m-1) = S(:,2:end);
    Y(:,end) = y;
    S(:,end) = s;
  end
  theta = (transpose(y)*y)/(transpose(y)*s);
  W = [Y theta*S];
  A = transpose(S)*Y;
  L = tril(A,-1);
  D = -1*diag(diag(A));
  MM = [D transpose(L); L theta*transpose(S)*S];
  M = inv(MM);
  
  % update the iteration
  k = k+1;
  if (xhistory)
    xhist = [xhist x];
  end
  
  % print some useful information
  if (display)
    opt = get_optimality(x,g,l,u);
    fprintf('%3d %16.8f %16.8f\n',k,f,opt);
  end
  
end

if (k == max_iters)
  fprintf(' warning: maximum number of iterations reached\n')
end

if ( get_optimality(x,g,l,u) < tol )
  fprintf(' stopping because convergence tolerance met!\n')
end

end
