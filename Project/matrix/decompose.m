function [decomposed_matrix] = decompose(matrix)
    decomposed_matrix = cell(size(matrix,1), size(matrix,2));
    for i=1:size(matrix,1)
        for j=1:size(matrix,2)
            decomposed_matrix{i,j} = squeeze(matrix(i,j,:));
        end
    end
end