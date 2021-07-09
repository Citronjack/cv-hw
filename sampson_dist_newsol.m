function sd = sampson_dist_newsol(F, x1_pixel, x2_pixel)
    % This function calculates the Sampson distance based on the fundamental matrix F
        % This function calculates the Sampson distance based on the fundamental matrix F
    e3 = [0,-1,0;1,0,0;0,0,0];
    sd = zeros([1, size(x1_pixel,2)]);
    i = 1;
    

    nom = (x2_pixel(:,i)'*F*x1_pixel(:,i)).^2;
    denom = (norm(e3*F*x1_pixel(:,i))).^2 + (norm(x2_pixel(:,i)'*F*e3)).^2;
    sd(i) = nom/denom;
    i = i+1;
    if i <= size(x1_pixel,2)
        GO TO 9
    end
end