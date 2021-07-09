function [correspondences_robust, largest_set_F] = F_ransac(correspondences, varargin)
    % This function implements the RANSAC algorithm to determine 
    % robust corresponding image points
   
    t = inputParser;
    valid_epsilon = @(x) isnumeric(x)  && (x<1) && (x>0);
    valid_p = @(x) isnumeric(x) && (x<1) && (x>0);
    valid_tolerance = @(x) isnumeric(x);

    addRequired(t, 'correspondences')
    addOptional(t,'epsilon', 0.5, valid_epsilon);
    addOptional(t,'p', 0.5, valid_p);
    addOptional(t,'tolerance', 0.01, valid_tolerance);
    
    parse(t, correspondences, varargin{:});
    epsilon = t.Results.epsilon;
    p =  t.Results.p;
    tolerance = t.Results.tolerance;  

    
    %% x1_pixel and x2_pixel from corresponcences and trafo into homogenous
    x1_pixel = ones([3, size(correspondences,2)]);
    x2_pixel = ones([3, size(correspondences,2)]);
    x1_pixel(1:2,:) = correspondences(1:2,:);
    x2_pixel(1:2,:) = correspondences(3:4,:);
    
    %% Preparation of cariables

    k = 8;
    s = log(1-p)/(log(1-(1-epsilon)^k));
    
    largest_set_size = 0;
    largest_set_dist = inf;
    largest_set_F = zeros([3,3]);
    
    %% RanSaC Algorithm

    for i = 1:k
       pos = randi(size(x1_pixel,2), 1,k);
       F = epa(correspondences(:, pos));
       
       dist = sampson_dist(F, x1_pixel, x2_pixel);
       subs = find(dist < tolerance);
       cons_set = [correspondences(:, subs); dist(subs)];
       
       set_cost = sum(cons_set(end, :));
       set_size = size(cons_set,2);
       
       if largest_set_size == 0
           largest_set_size = set_size;
           largest_set_dist = set_cost;
           best_cons_set = cons_set;
       elseif ~isempty(subs)
           if set_size > largest_set_size
               largest_set_size = set_size;
               largest_set_dist = set_cost;
               largest_set_F = F;
               best_cons_set = cons_set;
           elseif set_size == largest_set_size
               if set_cost < largest_set_dist
                    largest_set_size = set_size;
                    largest_set_dist = set_cost;
                    largest_set_F = F;
                    best_cons_set = cons_set;
               end
           end
       end
       
       
    end
    correspondences_robust = best_cons_set(1:4,:);
end