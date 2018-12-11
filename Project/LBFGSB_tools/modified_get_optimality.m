function [opt] = get_optimality(x,g,l,u,history_x,c,i,j)
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
   
   if(c < 2)
    opt = 9999;
   else
    ma = 0;
    for k=1:c-1
        ma = ma + abs(max(history_x{1,k+1}{i,j}-history_x{1,k}{i,j}));
    end
    opt = ma/c;
   end
  
end