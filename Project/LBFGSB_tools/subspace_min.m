function [xbar, line_search_flag] = subspace_min(x,g,l,u,xc,c,theta,W,M)
   % function [xbar] = subspace_min(x,g,l,u,xc,c,theta,W,M)
   % Subspace minimization for the quadratic model over free variables.
   % Direct Primal Method, Page 12.
   % INPUTS:
   %  x: [n,1] design vector.
   %  g: [n,1] objective gradient vector.
   %  l: [n,1] lower bound vector.
   %  u: [n,1] uppder bound vector.
   %  xc: [n,1] generalized Cauchy point.
   %  c: [2m,1] minimization initialization vector.
   %  theta: positive LBFGS scaling parameter.
   %  W: [n,2m] LBFGS matrix storage.
   %  M: [2m,2m] LBFGS matrix storage.
   
   % set the line search flag to true
   line_search_flag = true;
   
   % compute the free variables.
   n = length(x);
   free_vars_idx = [];
   Z = [];
   for i=1:length(xc)
     if ( ( xc(i) ~= u(i) ) && ( xc(i) ~= l(i) ) )
       free_vars_idx = [free_vars_idx i];
       unit = zeros(n,1);
       unit(i) = 1;
       Z = [Z unit];
     end
   end
   num_free_vars = length(free_vars_idx);
   
   if (num_free_vars == 0)
     xbar = xc;
     line_search_flag = false;
     return;
   end
   
   % compute W^T Z, the restriction of W to free variables.
   WTZ = transpose(W)*Z;
   
   % compute the reduced gradient of mk restricted to free variables.
   rr = g + theta*(xc-x) - W*(M*c);
   r = zeros(num_free_vars, 1);
   for i=1:num_free_vars
     r(i) = rr(free_vars_idx(i));
   end
   
   % form intermediate variables.
   invtheta = 1.0/theta;
   v = M*(WTZ*r);
   N = invtheta*WTZ*transpose(WTZ);
   N = eye(size(N)) - M*N;
   v = N\v;
   du = -invtheta*r - invtheta^2 * transpose(WTZ)*v;
   
   % find alpha star
   alpha_star = find_alpha(l,u,xc,du,free_vars_idx);
   
   % compute the subspace minimization
   d_star = alpha_star*du;
   xbar = xc;
   for i=1:num_free_vars
     idx = free_vars_idx(i);
     xbar(idx) = xbar(idx) + d_star(i);
   end
   
end