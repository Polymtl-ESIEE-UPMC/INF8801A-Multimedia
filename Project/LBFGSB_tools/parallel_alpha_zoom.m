function [alpha] = parallel_alpha_zoom(func,x0,f0,g0,p,alpha_lo,alpha_hi,i,j,alpha)
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
    k = 0;
    max_iters = 20;
    dphi0 = transpose(g0{i,j})*p{i,j};
    x = x0;
    
     while true
       alpha_i = 0.5*(alpha_lo{i,j} + alpha_hi{i,j});
       alpha{i,j} = alpha_i;
       x{i,j} = x0{i,j} + alpha_i*p{i,j};
       [f_i,g_i] = feval(func,x);
       x_lo = x0;
       x_lo{i,j} = x0{i,j} + alpha_lo{i,j}*p{i,j};
       f_lo = feval(func, x_lo);
       
       if ( (norm(f_i{i,j},size(f_i{i,j},1)) > norm(f0{i,j},size(f_i{i,j},1)) + c1*alpha_i*dphi0) || ( norm(f_i{i,j},size(f_i{i,j},1)) >= norm(f_lo{i,j},size(f_i{i,j},1))) )
         alpha_hi{i,j} = alpha_i;
       else
         dphi = transpose(g_i{i,j})*p{i,j};
         if ( ( abs(dphi) <= -c2*dphi0 ) )
           alpha{i,j} = alpha_i;
           break;
         end
         if ( dphi * (alpha_hi{i,j}-alpha_lo{i,j}) >= 0 )
           alpha_hi{i,j} = alpha_lo{i,j};
         end
         alpha_lo{i,j} = alpha_i;
       end
       k = k+1;
       if (k > max_iters)
         alpha{i,j} = alpha_i;
         break;
       end
     end
 end