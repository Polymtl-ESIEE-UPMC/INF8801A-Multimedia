function [x] = wrapper_strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached)
   [alpha] = strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached);
   for i=1:size(x,1)
      for j=1:size(x,2)
         x{i,j} = x{i,j} + alpha{i,j} * (xbar{i,j} - x{i,j});
      end
   end
end