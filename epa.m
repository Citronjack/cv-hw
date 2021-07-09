function EF = epa(correspondences, K)
    % Depending on whether a calibrating matrix 'K' is given,
    % this function calculates either the essential or the fundamental matrix
    % with the eight-point algorithm.
 
    
    %% First part of the Code
    x1= ones(3, size(correspondences, 2));
    x2= ones(3, size(correspondences, 2));
    
    x1(1:2, :) = correspondences(1:2, :);
    x2(1:2, :) = correspondences(3:4, :);
    x1 = x1';
    x2 = x2';

    % Camera calibration
    if exist('K', 'var')
        K = inv(K);
        for i = 1:size(x1,1)
            x1(i,:) = transpose(K*x1(i,:)');
            x2(i,:) = transpose(K*x2(i,:)');
        end
    end
    
    
    A = zeros([size(x1,1), 9]);
    % Kronecker Produkt to calculate A
    for i = 1:size(A,1)
       A(i,:) =  kron(x1(i,:), x2(i,:));
    end
    x1 = x1';
    x2 = x2';
    % SVD to calculate orthogonal matrix V
    [U,S,V] = svd(A);

    %% Estimation of the matrices
    % Take last column and trafo into 3x3
    %% Estimation of the matrices
    % Take last column and trafo into 3x3
    
    G = reshape(V(:, end), [3,3]);
    % find essential matrix E
    [G_U, G_S, G_V] = svd(G);
    T = eye(3);
    T(end,end) = 0;
    if exist('K', 'var')
        EF = G_U*T*G_V';
    else
        G_S(end,end) = 0;
        EF = G_U*G_S*G_V';
    end

    
end

