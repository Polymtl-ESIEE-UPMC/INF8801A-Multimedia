function [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached,history_x,c)
   opt = cell(size(x,1), size(x,2));
   all_optimal_reached = true;
   for i=1:size(x,1)
      for j=1:size(x,2)
         opt{i,j} = modified_get_optimality(x{i,j},g{i,j},l,u,history_x,c,i,j);
         if( optimal_reached{i,j} == false )
            if (opt{i,j} < tol)
               optimal_reached{i,j} = true;
            else
               all_optimal_reached = false;   
            end
         end
      end
   end
end