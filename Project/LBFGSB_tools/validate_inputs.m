function [x0] = validate_inputs(func,x0,l,u)
   % function [x0] = validate_inputs(func,x0,l,u)
   % Validate the inputs to the LBFGSB algorithm.
   % INPUTS:
   %  func: function handle to the objective, returns [f,g].
   %  x0: [n,1] initial design vector.
   %  l: [n,1] lower bound vector.
   %  u: [n,1] upper bound vector.
   % OUTPUTS:
   %  none
   
   % sanity check of inputs
   if ( nargout(func) ~=2 && nargout(func) ~= -1 )
     error('input func return must be of form [f,g]');
   end
   if ( ~ iscolumn(x0) )
     error('input x0 must be a column vector')
   end
   if ( ~ iscolumn(l) )
     error('input l must be a column vector')
   end
   if ( ~ iscolumn(u) )
     error('input u must be a column vector')
   end
   if ( length(l) ~= length(x0) )
     error('input l must be of equal length to input x0');
   end
   if ( length(u) ~= length(x0) )
     error('input u must be of equal length to input xo');
   end
   
   % pull back x0 into the feasible design space if needed
   modified_x0 = false;
   for i=1:length(x0)
     if ( x0(i) < l(i) )
       x0(i) = l(i);
       modified_x0 = true;
     elseif ( x0(i) > u(i) )
       x0(i) = u(i);
       modified_x0 = true;
     end
   end
   if (modified_x0)
     fprintf(' note: initial guess x0 outside of bounds\n');
     fprintf('       projecting x0 back into the feasible space\n');
   end
   
end