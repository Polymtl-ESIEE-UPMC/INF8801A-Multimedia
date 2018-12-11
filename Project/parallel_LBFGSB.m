function [x,xhist] = parallel_LBFGSB(func,x0,l,u,options)
    warning('off')
   % function [x,xhist] = LBFGSB(func,x0,l,u,options)
   % Perform bound-constrained optimization with L-BFGS-B.
   % INPUTS:
   %  x0: [n,1] initial design vector.
   %  l: [n,1] lower bound constraint vector.
   %  u: [n,1] upper bound constraint vector.
   %  options: matlab struct with the optional components:
   %   'm': the maximum number of stored L-BFGS iteration pairs.
   %   'tol': the convergence tolerance for the projected gradient.
   %   'display': true/false - should iterations be displayed?
   %   'xhist': true/false - should the entire search history be stored?
   
   
   % validate the inputs
   [x0] = wrapper_validate_inputs(func,x0,l,u);
   
   % set options
   [m,tol,max_iters,display,xhistory] = set_options(options);
   
   % initialize BFGS variables
   n = length(x0{1,1});

   Y = cell(size(x0,1), size(x0,2));
   S = cell(size(x0,1), size(x0,2));
   W = cell(size(x0,1), size(x0,2));
   M = cell(size(x0,1), size(x0,2));
   theta = cell(size(x0,1), size(x0,2));
   optimal_reached = cell(size(x0,1), size(x0,2));
   alpha = cell(size(x0,1), size(x0,2));
   for i=1:size(x0,1)
    for j=1:size(x0,2)
      Y{i,j} = zeros(n,0);
      S{i,j} = zeros(n,0);
      W{i,j} = zeros(n,1);
      M{i,j} = zeros(1,1);
      theta{i,j} = 1;
      optimal_reached{i,j} = false;
      alpha{i,j} = max_iters;
    end
   end
   history_x = cell(1,10);
   for i=1:size(history_x,2)
    history_x{1,i} = cell(size(x0,1), size(x0,2));
   end
   history_x{1,1} = zeros(n,1);
   check_index = 0;
   

   % initialize objective variables
   x = x0;
  
   [f,g] = feval(func, x);
   
   % initialize quasi-Newton iterations
   k = 0;
   
   % print out some useful information, if specified
   if (display)
     fprintf(' iter        f(x)          optimality\n')
     fprintf('-------------------------------------\n')
     opt = get_optimality(x,g,l,u);
     fprintf('%3d %16.8f %16.8f\n',k,f,opt);
   end
   
   % save the xhistory, if specified
   xhist = [];
   if (xhistory)
     xhist = [xhist x0];
   end

   [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached,history_x,check_index);
   % perform quasi-Newton iterations
   while ( (all_optimal_reached == false) && (k < max_iters) )
     fprintf(' iteration: %d\n',k);
     fprintf('>=');
     
     % update search information
     x_old = x;
     g_old = g;
     
     % compute the new search direction
     [xc, c] = wrapper_get_cauchy_point(x,g,l,u,theta{i,j},W{i,j},M{i,j},optimal_reached);
     fprintf('=');
     [xbar, line_search_flag] = wrapper_subspace_min(x,g,l,u,xc,c,theta{i,j},W{i,j},M{i,j},optimal_reached);
     fprintf('=');
     [x] = wrapper_strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached,alpha);
     fprintf('=');

     % update the LBFGS data structures
     [f,g] = feval(func, x);
     
     y = cell(size(x0,1), size(x0,2));
     s = cell(size(x0,1), size(x0,2));
     for i=1:size(x0,1)
      for j=1:size(x0,2)
        y{i,j} = g{i,j} - g_old{i,j};
        s{i,j} = x{i,j} - x_old{i,j};
      end
     end     

     curv = cell(size(x0,1), size(x0,2));
     for i=1:size(x0,1)
      for j=1:size(x0,2)
        curv{i,j} = abs(transpose(s{i,j})*y{i,j});
        if (curv{i,j} < eps)
          fprintf(' warning: negative curvature detected\n');
          fprintf('          skipping L-BFGS update\n');
          %k = k+1;
          %continue;
        else
          if (k < m)
            Y{i,j} = [Y{i,j} y{i,j}];
            S{i,j} = [S{i,j} s{i,j}];
          else
            size(Y{i,j})
            size(Y{i,j}(:,1:m-1))
            size(Y{i,j}(:,2:end))
            if(size(Y{i,j},2) >= m)
                Y{i,j}(:,1:m-1) = Y{i,j}(:,2:end);
                S{i,j}(:,1:m-1) = S{i,j}(:,2:end);
                Y{i,j}(:,end) = y{i,j};
                S{i,j}(:,end) = s{i,j};
            else
                Y{i,j}(:,1:end-1) = Y{i,j}(:,2:end);
                S{i,j}(:,1:end-1) = S{i,j}(:,2:end);
                Y{i,j}(:,end) = y{i,j};
                S{i,j}(:,end) = s{i,j};
            end

          end
          theta{i,j} = (transpose(y{i,j})*y{i,j})/(transpose(y{i,j})*s{i,j});
          W{i,j} = [Y{i,j} theta{i,j}*S{i,j}];
          A = transpose(S{i,j})*Y{i,j};
          L = tril(A,-1);
          D = -1*diag(diag(A));
          MM = [D transpose(L); L theta{i,j}*transpose(S{i,j})*S{i,j}];
          M{i,j} = inv(MM);
        end
      end
     end
     fprintf('=');
     
     % update the iteration
     k = k+1;
     if (xhistory)
       xhist = [xhist x];
     end
     
     if(k>1)
         if(check_index == 10)
            check_index = 0;
            for i=1:check_index-1
                history_x{1,i} = history_x{1,i+1}
            end
         else
            check_index = check_index+1; 
         end
         history_x{1,check_index} = x;
     end
     
     [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached,history_x,check_index);
     
     % print some useful information
     if (display)
       [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached);
       fprintf('%3d %16.8f %16.8f\n',k,f,opt);
     end
     fprintf('=<\n');

   end
   
   if (k == max_iters)
     fprintf(' warning: maximum number of iterations reached\n')
   end
   
   if ( all_optimal_reached == true )
     fprintf(' stopping because convergence tolerance met!\n')
   end
   
end
   