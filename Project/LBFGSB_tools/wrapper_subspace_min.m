function [xbar, line_search_flag] = wrapper_subspace_min(x,g,l,u,xc,c,theta,W,M,optimal_reached)
   xbar = cell(size(x,1), size(x,2));
   line_search_flag = cell(size(x,1), size(x,2));
   for i=1:size(x,1)
      for j=1:size(x,2)
         if(optimal_reached{i,j} == false)
            [xbar{i,j}, line_search_flag{i,j}] = subspace_min(x{i,j},g{i,j},l,u,xc{i,j},c{i,j},theta,W,M);
         end
      end
   end
end