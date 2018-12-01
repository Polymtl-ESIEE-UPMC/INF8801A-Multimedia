function [t,d,F] = get_breakpoints(x,g,l,u)
   % function [t,d,F] = get_breakpoints(x,g,l,u)
   % Compute the breakpoint variables needed for the Cauchy point.
   % Equations (4.1),(4.2), and F in Algorigthm CP: Initialize.
   % INPUTS:
   %  x: [n,1] design vector.
   %  g: [n,1] objective gradient.
   %  l: [n,1] lower bound vector.
   %  u: [n,1] upper bound vector.
   % OUTPUTS:
   %  t: [n,1] breakpoint vector.
   %  d: [n,1] search direction vector.
   %  F: [n,1] the indices that sort t from low to high.
   
   n = length(x);
   t = zeros(n,1);
   d = -g;
   for i=1:n
     if ( g(i) < 0 )
       t(i) = ( x(i) - u(i) ) / g(i);
     elseif ( g(i) > 0 )
       t(i) = ( x(i) - l(i) ) / g(i);
     else
       t(i) = realmax;
     end
     if ( t(i) < eps )
       d(i) = 0.0;
     end
   end
   tuple = [t linspace(1,n,n)'];
   tuple = sortrows(tuple);
   F = tuple(:,2);
   
end