function [matrix] = compose(decomposed_matrix)
    matrix = zeros(size(decomposed_matrix, 1), size(decomposed_matrix, 2), size(decomposed_matrix{1,1},1));
    for d=1:size(decomposed_matrix{1,1},1)
        for i=1:size(decomposed_matrix, 1)
            for j=1:size(decomposed_matrix, 2)
                matrix(i,j,d) = decomposed_matrix{i,j}(d);
            end
        end
    end
end