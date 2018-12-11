function [x] = wrapper_strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached,alpha)
%    [alpha] = parallel_strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached);
   for i=1:size(x,1)
      for j=1:size(x,2)
%          for m=1:size(xbar{i,j},1)
%              for n=1:size(xbar{i,j},2)
%                 if(xbar{i,j}(m,n) == -Inf || xbar{i,j}(m,n) == Inf)
%                     xbar{i,j}(m,n) = 0;
%                 end
%              end
%          end
         local_alpha = 0;
         if(optimal_reached{i,j} == false && alpha{i,j} ~= 0)
            if(line_search_flag{i,j} == false)
                local_alpha = 1;
            else
                local_alpha = alpha{i,j};
                alpha{i,j} = alpha{i,j}-1;
            end
         end
         inverse = -1;
         if(local_alpha == 0)
             inverse = 1;
         end
         x{i,j} = (x{i,j} - local_alpha * (xbar{i,j} - x{i,j}));
      end
   end
end