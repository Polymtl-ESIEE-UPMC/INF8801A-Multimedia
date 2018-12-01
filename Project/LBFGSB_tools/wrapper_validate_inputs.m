function [x] = wrapper_validate_inputs(func,x,l,u)
   for i=1:size(x,1)
      for j=1:size(x,2)
        x{i,j} = validate_inputs(func,x{i,j},l,u);
      end
   end
end