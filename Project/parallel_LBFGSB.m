function [x,xhist] = parallel_LBFGSB(func,x0,l,u,options)
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
   Y = zeros(n,0);
   S = zeros(n,0);
   W = zeros(n,1);
   M = zeros(1,1);
   theta = 1;
   
   % initialize objective variables
   x = x0;
   f = cell(size(x0,1), size(x0,2));
   g = cell(size(x0,1), size(x0,2));
   
   % TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    for i=1:size(x0,1)
    %      for j=1:size(x0,2)
    %        [f{i,j},g{i,j}] = feval(func, x{i,j});
    %      end
    %    end
        [f,g] = feval(func, x);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
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

   optimal_reached = cell(size(x,1), size(x,2));
   for i=1:size(optimal_reached,1)
    for j=1:size(optimal_reached,2)
      optimal_reached{i,j} = false;
    end
   end
   [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached);
   % perform quasi-Newton iterations
   while ( (all_optimal_reached == false) && (k < max_iters) )
     
      
     % update search information
     x_old = x;
     g_old = g;
     
     % compute the new search direction
     [xc, c] = wrapper_get_cauchy_point(x,g,l,u,theta,W,M,optimal_reached);
     [xbar, line_search_flag] = wrapper_subspace_min(x,g,l,u,xc,c,theta,W,M,optimal_reached);
   
     [x] = wrapper_strong_wolfe(func,x,f,g,xbar,line_search_flag,optimal_reached);
     
     % update the LBFGS data structures
     % TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    for i=1:size(x0,1)
%      for j=1:size(x0,2)
%        [f{i,j},g{i,j}] = feval(func, x{i,j});
%      end
%    end
    [f,g] = feval(func, x);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      y = cell(size(x0,1), size(x0,2));
      s = cell(size(x0,1), size(x0,2));
     for i=1:size(x0,1)
      for j=1:size(x0,2)
        y{i,j} = g{i,j} - g_old{i,j};
        s{i,j} = x{i,j} - x_old{i,j};
      end
     end     

     % TODO: HOW TO CHANGE PARAMETER ? %%%%%%%%%%%
     curv = abs(transpose(s)*y);
     if (curv < eps)
       fprintf(' warning: negative curvature detected\n');
       fprintf('          skipping L-BFGS update\n');
       k = k+1;
       continue;
     end
     if (k < m)
       Y = [Y y];
       S = [S s];
     else
       Y(:,1:m-1) = Y(:,2:end);
       S(:,1:m-1) = S(:,2:end);
       Y(:,end) = y;
       S(:,end) = s;
     end
     theta = (transpose(y)*y)/(transpose(y)*s);
     W = [Y theta*S];
     A = transpose(S)*Y;
     L = tril(A,-1);
     D = -1*diag(diag(A));
     MM = [D transpose(L); L theta*transpose(S)*S];
     M = inv(MM);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     % update the iteration
     k = k+1;
     if (xhistory)
       xhist = [xhist x];
     end
     
     % print some useful information
     if (display)
       [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached);
       fprintf('%3d %16.8f %16.8f\n',k,f,opt);
     end
     
   end
   
   if (k == max_iters)
     fprintf(' warning: maximum number of iterations reached\n')
   end
   
   [opt, optimal_reached, all_optimal_reached] = wrapper_get_optimality(x,g,l,u,tol,optimal_reached);
   if ( all_optimal_reached == true )
     fprintf(' stopping because convergence tolerance met!\n')
   end
   
   end
   