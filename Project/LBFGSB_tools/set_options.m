function [m,tol,max_iters,display,xhistory] = set_options(options)
   % function [m,tol,max_iters] = set_options(options)
   % Set optionally defined user input parameters.
   % INPUTS:
   %  options: a matlab struct ('m','tol','max_iters').
   % OUTPUTS:
   %  m: the limited memory storage size.
   %  tol: the convergence tolerance criteria.
   %  max_iters: the maximum number of quasi-Newton iterations.
   %  display: true/false should iteration information be displayed?
   %  xhistory: true/false should the entire search history be stored?
   m = 10;
   tol = 1.0e-5;
   max_iters = 20;
   display = false;
   xhistory = false;
   if ( isfield(options, 'm') )
     m = options.m;
   end
   if ( isfield(options, 'tol') )
     tol = options.tol;
   end
   if ( isfield(options, 'max_iters') )
     max_iters = options.max_iters;
   end
   if ( isfield(options, 'display') )
     display = options.display;
   end
   if ( isfield(options, 'xhistory') )
     xhistory = options.xhistory;
   end
   
end