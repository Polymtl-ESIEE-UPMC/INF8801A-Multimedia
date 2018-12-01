function [alpha] = alpha_zoom(func,x0,f0,g0,p,alpha_lo,alpha_hi)
   % function [alpha] = alpha_zoom(func,x0,f0,g0,p,alpha_lo,alpha_hi)
   % Algorithm 3.6, Page 61. "Numerical Optimization". Nocedal & Wright.
   % INPUTS:
   %  func: objective function handle.
   %  x0: [n,1] initial design vector.
   %  f0: initial objective value.
   %  g0: [n,1] initial objective gradient vector.
   %  p: [n,1] search direction vector.
   %  alpha_lo: low water mark for alpha.
   %  alpha_hi: high water mark for alpha.
   % OUTPUTS:
   %  alpha: zoomed in alpha.
   
   % initialize variables
   c1 = 1e-4;
   c2 = 0.9;
   i = 0;
   max_iters = 20;
   dphi0 = transpose(g0)*p;
   
   while true
     alpha_i = 0.5*(alpha_lo + alpha_hi);
     alpha = alpha_i;
     x = x0 + alpha_i*p;
     [f_i,g_i] = feval(func,x);
     x_lo = x0 + alpha_lo*p;
     f_lo = feval(func, x_lo);
     if ( (f_i > f0 + c1*alpha_i*dphi0) || ( f_i >= f_lo) )
       alpha_hi = alpha_i;
     else
       dphi = transpose(g_i)*p;
       if ( ( abs(dphi) <= -c2*dphi0 ) )
         alpha = alpha_i;
         break;
       end
       if ( dphi * (alpha_hi-alpha_lo) >= 0 )
         alpha_hi = alpha_lo;
       end
       alpha_lo = alpha_i;
     end
     i = i+1;
     if (i > max_iters)
       alpha = alpha_i;
       break;
     end
   end
   
end