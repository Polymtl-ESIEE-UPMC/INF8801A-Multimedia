function [output] = wrapper_LBFGSB(input)

  global composed
  composed = input;
  decomposed = decompose(input);
  l = zeros(size(decomposed{1,1}));
  u = 1000.0*ones(size(decomposed{1,1}));
  
  x = cell(size(decomposed,1), size(decomposed,2));
  for i=1:size(x,1)
    for j=1:size(x,2)
      x{i,j} = ones(size(decomposed{i,j}));
    end
  end

%   [f,g] =lost_func_testing(x);
  x = parallel_LBFGSB(@lost_func_testing,x,l,u,[]);
  output = compose(x);
end

function [f,g] = lost_func(x)
  f = transpose(x)*x;
  if (nargout > 1)
    g = 2.0*x;
  end
end

function [f,g] = lost_func_testing(x)
  global composed
  composedx = compose(x);
  temp = composedx - composed;
  f = decompose(temp);
  if (nargout > 1)
    g = cell(size(f,1), size(f,2));
    for i=1:size(g,1)
      for j=1:size(g,2)
        g{i,j} = gradient(f{i,j});
      end
    end
  end
end