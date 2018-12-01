function [xc, c] = get_cauchy_point(x,g,l,u,theta,W,M)
   % function [xc, c] = get_cauchy_point(x,g,l,u,theta,W,M)
   % Compute the generalized Cauchy point.
   % Algorithm CP, Pages 8-9.
   % INPUTS:
   %  x: [n,1] design vector.
   %  g: [n,1] objective function gradient.
   %  l: [n,1] lower bound vector.
   %  u: [n,1] uppder bound vector.
   %  theta: positive BFGS scaling.
   %  W: [n,2m] BFGS matrix storage.
   %  M: [2m,2m] BFGSB matrix storage.
   % OUTPUTS:
   %  xc - [n,1] the generalized Cauchy point.
   %  c - [2m,1] initialization vector for subspace minimization.
   
   % perform the initialization step
   [tt,d,F] = get_breakpoints(x,g,l,u);
   xc = x;
   p = transpose(W) * d;
   c = zeros(size(W,2),1);
   fp = -transpose(d)*d;
   fpp = -theta*fp - transpose(p)*M*p;
   fpp0 = -theta*fp;
   dt_min = -fp/fpp;
   t_old = 0;
   for j=1:length(x)
     i = j;
     if (F(i) > 0)
       break;
     end
   end
   b = F(i);
   t = tt(b);
   dt = t-t_old;
   
   % examine the subsequent segments
   while ( (dt_min > dt) && (i <= length(x)) )
     if ( d(b) > 0)
       xc(b) = u(b);
     elseif ( d(b) < 0)
       xc(b) = l(b);
     end
     zb = xc(b) - x(b);
     c = c + dt*p;
     gb = g(b);
     wbt = W(b,:);
     fp = fp + dt*fpp + gb*gb + theta*gb*zb - gb*wbt*(M*c);
     fpp = fpp - theta*gb*gb - 2.0*gb*wbt*(M*p) - gb*gb*wbt*(M*transpose(wbt));
     fpp = max( eps*fpp0, fpp );
     p = p + gb*transpose(wbt);
     d(b) = 0.0;
     dt_min = -fp/fpp;
     t_old = t;
     i = i+1;
     if (i <= length(x))
       b = F(i);
       t = tt(b);
       dt = t - t_old;
     end
   end
   
   % perform final updates
   dt_min = max(dt_min, 0);
   t_old = t_old + dt_min;
   for j=i:length(xc)
     idx = F(j);
     xc(idx) = x(idx) + t_old*d(idx);
   end
   c = c + dt_min*p;
   
end