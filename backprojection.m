function [repro_error, x2_repro] = backprojection(correspondences, P1, Image2, T, R, K)
    % This function calculates the mean error of the back projection
    % of the world coordinates P1 from image 1 in camera frame 2
    % and visualizes the correct feature coordinates as well as the back projected ones.
    
    %world coordinates of camerframe 2
   
    P2 = R*P1 + T;

    
    x1 = ones(3, size(correspondences,2));
    x2 = ones(3, size(correspondences,2));
    x1(1:2,:) = correspondences(1:2,:);
    x2(1:2,:) = correspondences(3:4,:);
    %x1_k = K\x1;
    x2_rec = K*P2;
    lambda = x2_rec(3,:);
    x2_repro = x2_rec./lambda;
    
    imshow(Image2);hold on
    scatter(x2(1,:), x2(2,:), 'MarkerEdgeColor', 'c')
    text(x2(1,:), x2(2,:),num2str([1:size(x2,2)]')) 
    scatter(x2_repro(1,:), x2_repro(2,:), 'MarkerEdgeColor', 'r')
    text(x2_repro(1,:), x2_repro(2,:),num2str([1:size(x2,2)]')) 
    
    for i = 1:size(x2,2)
       plot([x2(1,i), x2_repro(1,i)], [x2(2,i), x2_repro(2,i)], 'g')
    end
    
    repro_error = norm(x2(1:3,:)-x2_repro(1:3,:),2)*1/(size(x2,2));
end

