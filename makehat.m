function [ cp_cell ] = makehat(mat_3xN)
% MAKEHAT generates a cell array whjich contains all the corss product
% matrices from an 3xN input matrix
    cp_cell = cell(1, size(mat_3xN,2));
    for i = 1:size(mat_3xN,2)
        vec_hat = zeros(3,3);
        vec_hat(1,2) = -mat_3xN(3, i);
        vec_hat(1,3) = mat_3xN(2, i);
        vec_hat(2,1) = mat_3xN(3, i);
        vec_hat(2,3) = -mat_3xN(1, i);
        vec_hat(3,1) = -mat_3xN(2, i);
        vec_hat(3,2) = mat_3xN(2, i);
        cp_cell{1,i} = vec_hat;
    end
end

