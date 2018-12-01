function [opt] = get_optimality(x,g,l,u)
   % function [opt] = get_optimality(x,g,l,u)
   % Get the inf-norm of the projected gradient.
   % Equation (6.1), Page 17.
   % INPUTS:
   %  x: [n,1] design vector.
   %  g: [n,1] objective function gradient.
   %  l: [n,1] lower bound vector.
   %  u: [n,1] upper bound vector.
   % OUTPUTS:
   %  opt: the inf-norm of the projected gradient.
   
   projected_g = x-g;
   for i=1:length(x)
     if (projected_g(i) < l(i))
       projected_g(i) = l(i);
     elseif (projected_g(i) > u(i))
       projected_g(i) = u(i);
     end
   end
   projected_g = projected_g - x;
   opt = max(abs(projected_g));
   
end