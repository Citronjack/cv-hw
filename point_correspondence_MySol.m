Image1 = imread('sceneL.png');
IGray1 = rgb_to_gray(Image1);

Image2 = imread('sceneR.png');
IGray2 = rgb_to_gray(Image2);
Image1 = imread('sceneL.png');
I1 = rgb_to_gray(Image1);

Image2 = imread('sceneR.png');
I2 = rgb_to_gray(Image2);
min_corr= 0.95;
window_length = 25;
do_plot = true;

Ftp1 = harris_detector(IGray1,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);
Ftp2 = harris_detector(IGray2,'segment_length',9,'k',0.05,'min_dist',50,'N',20,'do_plot',false);

Im1 = double(I1);
Im2 = double(I2);

    % In this function you are going to compare the extracted features of a stereo recording
    % with NCC to determine corresponding image points.
    
    
%     p = inputParser;
%     isodd = @(x) (mod(x,2) == 1);
%     valid_window_length = @(x) isnumeric(x) && isodd(x);
%     valid_min_corr = @(x) isnumeric(x) && (x<1) && (x>0);
%     valid_do_plot = @(x) islogical(x);
% 
%     addRequired(p, 'I1')
%     addRequired(p, 'I2')
%     addRequired(p, 'Ftp1')
%     addRequired(p, 'Ftp2')
%     addOptional(p,'window_length', 25, valid_window_length);
%     addOptional(p,'min_corr', 0.95, valid_min_corr);
%     addOptional(p,'do_plot', false, valid_do_plot);
%     
%     parse(p, I1, I2, Ftp1, Ftp2, varargin{:});
%     window_length = p.Results.window_length;
%     min_corr =  p.Results.min_corr;
%     do_plot = p.Results.do_plot;  
%     Im1 = double(I1);
%     Im2= double(I2);
    
    rows_I1 = size(I1, 1);
    cols_I1 = size(I1, 2);
    rows_I2 = size(I2, 1);
    cols_I2 = size(I2, 2);
    
    k1 = 1;
    for i = 1:size(Ftp1,2)
        % left side ok
        dist_left = Ftp1(1,i) - window_length/2;
        dist_right = Ftp1(1,i) + window_length/2;
        dist_up = Ftp1(2,i) - window_length/2;
        dist_bottom = Ftp1(2,i) + window_length/2;
        if (dist_left >= 0) && (dist_right <= cols_I1) && (dist_up >= 0) && (dist_bottom <= rows_I1)
            Ftp1_ok(k1) = i;  
            k1 = k1+1;
        end
    end
    
    % Checking if feature in I2 lie far enough away
    k2 = 1;
    for i = 1:size(Ftp2,2)
        % left side ok
        dist_left = Ftp2(1,i) - window_length/2;
        dist_right = Ftp2(1,i) + window_length/2;
        dist_up = Ftp2(2,i) - window_length/2;
        dist_bottom = Ftp2(2,i) + window_length/2;
        if (dist_left >= 0) && (dist_right <= cols_I1) && (dist_up >= 0) && (dist_bottom <= rows_I1)
            Ftp2_ok(k2) = i;  
            k2 = k2+1;
        end
    end
    Ftp1 = Ftp1(:,Ftp1_ok);
    Ftp2 = Ftp2(:,Ftp2_ok);
    no_pts1 = size(Ftp1,2);
    no_pts2 = size(Ftp2,2);
    
    
    dist = floor(window_length/2);
    k = 1;
    for i = 1:size(Ftp1, 2)
       x_coor =  Ftp1(1,i);
       y_coor =  Ftp1(2,i);
       x_start = x_coor - dist;
       x_end = x_coor + dist;
       
       y_start = y_coor - dist;
       y_end = y_coor + dist;
       
       W = Im1(y_start:y_end, x_start:x_end);
       
       W_std = std(W(:),0);
       W_mean = mean(W(:));
       W = 1/W_std*(W-W_mean);
       
       Mat_feat_1(:,k) = W(:);
       k = k+1;
    end
    
    k = 1;
    for i = 1:size(Ftp2, 2)
       x_coor =  Ftp2(1,i);
       y_coor =  Ftp2(2,i);
       x_start = x_coor - dist;
       x_end = x_coor + dist;
       
       y_start = y_coor - dist;
       y_end = y_coor + dist;
       
       W = Im2(y_start:y_end, x_start:x_end);
       
       W_std = std(W(:),0);
       W_mean = mean(W(:));
       W = 1/W_std*(W-W_mean);
       
       Mat_feat_2(:,k) = W(:);
       k = k+1;
    end
    
    
    % NCC calculatiuon old solution
    
%     NCC_matrix = zeros([size(Mat_feat_1,2), size(Mat_feat_2,2)]);
%     
%     for j = 1:size(Mat_feat_2,2)
%        im2_mat = reshape(Mat_feat_2(:,j), [window_length,window_length]);
% 
%         for i = 1:size(Mat_feat_1,2)
%             im1_mat = reshape(Mat_feat_1(:,i), [window_length,window_length]);
%             NCC_matrix(i,j) = 1/(window_length^2-1)*trace(im2_mat'*im1_mat);
%             
%         end
% 
%     end
%     %set all values lower than threshold to Zero!
%     NCC_matrix(NCC_matrix <= min_corr) = 0;
%     
%     %Sort in descending order
%     [NCC_sort, sorted_index] = sort(NCC_matrix(:), 'descend');
%     sorted_index = sorted_index(1:numel(NCC_sort(NCC_sort>0)));
   


    NCC_matrix = zeros([size(Mat_feat_1,2), size(Mat_feat_2,2)]);
    NCC_matrix = zeros([size(Mat_feat_2,2), size(Mat_feat_1,2)]);

    
    for j = 1:size(Mat_feat_1,2)
       im2_mat = reshape(Mat_feat_1(:,j), [window_length,window_length]);

        for i = 1:size(Mat_feat_2,2)
            im1_mat = reshape(Mat_feat_2(:,i), [window_length,window_length]);
            NCC_matrix(i,j) = 1/(window_length^2-1)*trace(im2_mat'*im1_mat);
            
        end

    end
    %set all values lower than threshold to Zero!
    NCC_matrix(NCC_matrix <= min_corr) = 0;
    [NCC_sort, sorted_index] = sort(NCC_matrix(:), 'descend');
    sorted_index = sorted_index(1:numel(NCC_sort(NCC_sort>0)));
    %Sort in descending order
    k = 1;
    for i= 1:length(sorted_index)
        ncc_idx = sorted_index(i);
        [row, col] = ind2sub(size(NCC_matrix), ncc_idx);

        if ~isnan(NCC_matrix(row, col))
           cor(:,k) = [Ftp1(:,col); Ftp2(:,row)];
           k = k+1;
           NCC_matrix(:, col) = NaN;
        end
        

    end
    
    
    %plot
    imshow(I1);
    hold on;
    transparent = imshow(I2);
    set(transparent, 'AlphaData', 0.5);
    scatter(cor(1,:), cor(2,:), 'MarkerEdgeColor', 'c')
    scatter(cor(3,:), cor(4,:), 'MarkerEdgeColor', 'r')
    
    for i = 1:size(cor,2)
       plot([cor(1,i), cor(3,i)], [cor(2,i), cor(4,i)], 'g')
    end


    
    
