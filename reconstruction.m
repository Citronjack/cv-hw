function [T, R, lambda, P1, camC1, camC2] = reconstruction(T1, T2, R1, R2, correspondences, K)
    %% Preparation
    % Init Cell arrays T_cell and R_cell
    T_cell = {T1, T2, T1, T2};
    R_cell = {R1, R1, R2, R2};
    
    % init homogenous coordinates x1 and x2
    x1 = ones(3, size(correspondences,2));
    x2 = ones(3, size(correspondences,2));
    x1(1:2,:) = correspondences(1:2,:);
    x2(1:2,:) = correspondences(3:4,:);
    x1 = K\x1;
    x2 = K\x2;
    
    %init d_cell
    d_cell = {zeros(size(correspondences,2),2), zeros(size(correspondences,2),2), zeros(size(correspondences,2),2), zeros(size(correspondences,2),2)};
    N = size(correspondences,2);
    positive_entries = zeros(1,3);
    %% Reconstruction
    % Constructing M1 and M2 as well as solving it via svd
    for i=1:4
        T = cell2mat(T_cell(i));
        R = cell2mat(R_cell(i));
        
        x1_hat_cell = makehat(x1);
        x2_hat_cell = makehat(x2);
        
        M1 = zeros(3*N, N+1);
        M2 = zeros(3*N, N+1);
        k = 1;
        for j = 2:3:3*N
           M1(j-1:j+1,k) = cell2mat(x2_hat_cell(k))*R*x1(:,k);
           M1(j-1:j+1, end) = cell2mat(x2_hat_cell(k))*T;
           
           M2(j-1:j+1,k) = cell2mat(x1_hat_cell(k))*R'*x2(:, k);
           M2(j-1:j+1, end) = -cell2mat(x1_hat_cell(k))*R'*T;
           k = k+1;
        end
        
        [~, ~, V_m1] = svd(M1);
        lambda_m1 = V_m1(:,end);
        lambda_m1  = lambda_m1/lambda_m1(end);
        [~, ~, V_m2] = svd(M2);
        lambda_m2 = V_m2(:,end);
        lambda_m2  = lambda_m2/lambda_m2(end);
        
        d_cell{i} = [lambda_m1(1:end-1), lambda_m2(1:end-1)];
        positive_entries(i) = length(find(lambda_m1(end-1)>0)) + length(find(lambda_m2(end-1)>0));
    end
    
    l = find(positive_entries == max(positive_entries));
    lambda = d_cell{l};
    
    T = cell2mat(T_cell(l));
    R = cell2mat(R_cell(l));
    
    %% Plotting
    P1 = lambda(:,1)'.*(x1);
    P2 = lambda(:,2)'.*(x2);
    
    camC1 = [[-0.2;0.2;1], [0.2;0.2;1], [0.2;-0.2;1], [-0.2;-0.2;1]];
    camC2 = R\(camC1 - T);
    
    %camPos1 = ;
    %camPos2 = R\(camPos1 - T);
    
    grid on
    hold on
    
    scatter3(P1(1,:), P1(2,:), P1(3,:),'k');
    text(P1(1,:), P1(2,:), P1(3,:),num2str([1:size(P1,2)]')); %axis([-0.5, 0.5, -0.8, 0.5, 0, 5]); 
    scatter3(P2(1,:), P2(2,:), P2(3,:),'k');
    text(P2(1,:), P2(2,:), P2(3,:),num2str([1:size(P2,2)]')); 
    
    
    %patch([-0.2,0.2,1,1], [0.2,0.2,1,1], [0.2,-0.2,1,1], [-0.2,-0.2,1,1],'FaceColor','r')
    campos([43;-22;-87]);
    camup([0 -1 0]);
    %patch(camC2(:,1), camC2(:,2), camC2(:,3), camC2(:,4),'FaceColor','r');
    
    plot3(camC1(1), camC1(2), camC1(3), 'Color', 'r', 'Marker', 's','MarkerSize', 20);
    text(camC1(1), camC1(2), camC1(3), 'Camera 1');
    plot3(camC2(1), camC2(2), camC2(3), 'Color', 'b', 'Marker', 's','MarkerSize', 20);
    text(camC2(1), camC2(2), camC2(3), 'Camera 2');
    
    xlabel('x')
    ylabel('x')
    zlabel('z')
    
    
end