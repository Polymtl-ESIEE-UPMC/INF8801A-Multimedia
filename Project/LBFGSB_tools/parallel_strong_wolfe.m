function [alpha] = parallel_strong_wolfe(func,x0,f0,g0,xbar,line_search_flag,optimal_reached)
   % function [alpha] = strong_wolfe(func,x0,f0,g0,p)
   % Compute a line search to satisfy the strong Wolfe conditions.
   % Algorithm 3.5. Page 60. "Numerical Optimization". Nocedal & Wright.
   % INPUTS:
   %  func: objective function handle.
   %  x0: [n,1] initial design vector.
   %  f0: initial function evaluation.
   %  g0: [n,1] initial objective gradient vector.
   %  p: [n,1] search direction vector.
   % OUTPUTS:
   % alpha: search length
   
   % initialize variables
%    c1 = 1e-4;
%    c2 = 0.9;
%    k = 0;
%    max_iters = 20;
%    alpha_max = 2.5;
% 
%    p = cell(size(x0,1), size(x0,2));
%    alpha = cell(size(x0,1), size(x0,2));
%    alpha_im1 = cell(size(x0,1), size(x0,2));
%    alpha_i = cell(size(x0,1), size(x0,2));
%    f_im1 = cell(size(x0,1), size(x0,2));
%    dphi0 = cell(size(x0,1), size(x0,2));
%    done = cell(size(x0,1), size(x0,2));
%    allDone = false;
% 
%    for i=1:size(x0,1)
%     for j=1:size(x0,2)
%       alpha{i,j} = 1;
%       p{i,j} = xbar{i,j} - x0{i,j};
%       alpha_im1{i,j} = 0;
%       alpha_i{i,j} = 1;
%       f_im1{i,j} = f0;
%       dphi0{i,j} = transpose(g0{i,j})*p{i,j}; 
%       done{i,j} = false;
%     end
%    end
%    
%    % search for alpha that satisfies strong-Wolfe conditions
%    while (allDone == false)
%      
%      x = cell(size(x0,1), size(x0,2));
%      for i=1:size(x0,1)
%       for j=1:size(x0,2)
%         x{i,j} = x0{i,j} + alpha_i{i,j}*p{i,j}; 
%       end
%      end
%      
%      [f_i,g_i] = feval(func, x);
% 
%      allDone = true;
% %      fprintf(' iteration: %d\n',k);
%      for i=1:size(x0,1)
%       for j=1:size(x0,2)
% %           fprintf(' i: %d ',i);
% %           fprintf(' j: %d\n',j);
%         if (done{i,j} == false && optimal_reached{i,j} == false && line_search_flag{i,j} == true)
%           if (norm(f_i{i,j},size(f_i{i,j},1)) > norm(f0{i,j},size(f_i{i,j},1)) + c1*dphi0{i,j}) || ( (k > 1) && (norm(f_i{i,j},size(f_i{i,j},1)) >= norm(f_im1{i,j},size(f_i{i,j},1))) )
%             alpha = parallel_alpha_zoom(func,x0,f0,g0,p,alpha_im1,alpha_i,i,j,alpha);
%             done{i,j} = true;
%             continue;
%           end
%           dphi = transpose(g_i{i,j})*p{i,j};
%           if ( abs(dphi) <= -c2*dphi0{i,j} )
%             alpha{i,j} = alpha_i{i,j};
%             done{i,j} = true;
%             continue;
%           end
%           if ( dphi >= 0 )
%             alpha = parallel_alpha_zoom(func,x0,f0,g0,p,alpha_i,alpha_im1,i,j,alpha);
%             done{i,j} = true;
%             continue;
%           end
%           allDone = false;
%           % update
%           alpha_i{i,j} = alpha_i{i,j} + 0.8*(alpha_max-alpha_i{i,j});
%         end
%       end
%      end
%      
%      % update
%      alpha_im1 = alpha_i;
%      f_im1 = f_i;
%      
%      if (k > max_iters)
%        alpha = alpha_i;
%        break;
%      end
%      
%      k = k+1;
%      
%    end
   
   alpha = cell(size(x0,1), size(x0,2));
   
   for i=1:size(x0,1)
    for j=1:size(x0,2)
      alpha{i,j} = 50;
    
    end
   end
   
end