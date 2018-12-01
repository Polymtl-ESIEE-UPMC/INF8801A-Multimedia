function [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached)
   opt = cell(size(x,1), size(x,2));
   all_optimal_reached = true;
   for i=1:size(x,1)
      for j=1:size(x,2)
         opt{i,j} = get_optimality(x{i,j},g{i,j},l,u);
         if( optimal_reached{i,j} ~= true && opt{i,j} < tol )
            optimal_reached{i,j} = true;
         else
            optimal_reached{i,j} = false;
            all_optimal_reached = false;
         end
      end
   end
end