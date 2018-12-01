function [alpha_star] = find_alpha(l,u,xc,du,free_vars_idx)
   % function [alpha_star] = find_alpha(l,u,xc,du,free_vars_idx)
   % Equation (5.8), Page 8.
   % INPUTS:
   %  l: [n,1] lower bound constraint vector.
   %  u: [n,1] upper bound constraint vector.
   %  xc: [n,1] generalized Cauchy point.
   %  du: [num_free_vars,1] solution of unconstrained minimization.
   % OUTPUTS:
   %  alpha_star: positive scaling parameter.
   
   alpha_star = 1;
   n = length(free_vars_idx);
   for i=1:n
     idx = free_vars_idx(i);
     if (du(i) > 0)
       alpha_star = min(alpha_star, ( u(idx)-xc(idx) )/du(i) );
     else
       alpha_star = min(alpha_star, ( l(idx)-xc(idx) )/du(i) );
     end
   end
   
end