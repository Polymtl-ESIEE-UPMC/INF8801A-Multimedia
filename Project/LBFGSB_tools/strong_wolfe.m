function [alpha] = strong_wolfe(func,x0,f0,g0,p)
   % function [alpha] = strong_wolfe(func,x0,f0,g0,p)
   % Compute a line search to satisfy the strong Wolfe conditions.
   % Algorithm 3.5. Page 60. "Numerical Optimization". Nocedal & Wright.
   % INPUTS:
   %  func: objective function handle.
   %  x0: [n,1] initial design vector.
   %  f0: initial function evaluation.
   %  g0: [n,1] initial objective gradient vector.
   %  p: [n,1] search direction vector.
   % OUTPUTS:
   % alpha: search length
   
   % initialize variables
   c1 = 1e-4;
   c2 = 0.9;
   alpha_max = 2.5;
   alpha_im1 = 0;
   alpha_i = 1;
   f_im1 = f0;
   dphi0 = transpose(g0)*p;
   i = 0;
   max_iters = 20;
   
   % search for alpha that satisfies strong-Wolfe conditions
   while true
     
     x = x0 + alpha_i*p;
     [f_i,g_i] = feval(func, x);
     if (f_i > f0 + c1*dphi0) || ( (i > 1) && (f_i >= f_im1) )
       alpha = alpha_zoom(func,x0,f0,g0,p,alpha_im1,alpha_i);
       break;
     end
     dphi = transpose(g_i)*p;
     if ( abs(dphi) <= -c2*dphi0 )
       alpha = alpha_i;
       break;
     end
     if ( dphi >= 0 )
       alpha = alpha_zoom(func,x0,f0,g0,p,alpha_i,alpha_im1);
       break;
     end
     
     % update
     alpha_im1 = alpha_i;
     f_im1 = f_i;
     alpha_i = alpha_i + 0.8*(alpha_max-alpha_i);
     
     if (i > max_iters)
       alpha = alpha_i;
       break;
     end
     
     i = i+1;
     
   end
   
end