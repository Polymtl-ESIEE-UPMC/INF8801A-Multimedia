function [output] = wrapper_LBFGSB(input)
  load('FBp3mapped.mat');
  input = FBpWarpedtemp;
  
  global composed
  composed = input;
  decomposed = decompose(input);
  l = zeros(size(decomposed{1,1}));
  u = 100000.0*ones(size(decomposed{1,1}));
  
  a = 0;
  b = 100000;
  x = cell(size(decomposed,1), size(decomposed,2));
  for i=1:size(x,1)
    for j=1:size(x,2)
%       x{i,j} = (b-a).*rand(size(decomposed{i,j}))+a;
       x{i,j} = ones(size(decomposed{i,j}));
    end
  end

  x = parallel_LBFGSB(@lost_func,x,l,u,[]);
  output = compose(x);
end

function [f,g] = lost_func(x)

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