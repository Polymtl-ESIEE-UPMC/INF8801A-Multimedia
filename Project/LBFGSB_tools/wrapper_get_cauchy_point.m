function [xc, c] = wrapper_get_cauchy_point(x,g,l,u,theta,W,M,optimal_reached)
   xc = cell(size(x,1), size(x,2));
   c = cell(size(x,1), size(x,2));
   for i=1:size(x,1)
      for j=1:size(x,2)
         if(optimal_reached{i,j} == false)
            [xc{i,j}, c{i,j}] = get_cauchy_point(x{i,j},g{i,j},l,u,theta,W,M);
         end
      end
   end
end