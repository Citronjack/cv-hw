function [T1, R1, T2, R2, U, V] = TR_from_E(E)
    % This function calculates the possible values for T and R 
    % from the essential matrix
    %% SVD von der essentiellen Matrix
    [U, S, V] = svd(E);
    
    %% U und V auf Rotationsmatrix überprüfen
    I = eye(3);
    I(3,3) = -1;
    if (det(U) ~= 1)
        disp('U is not a rotations matrix')
        U = U*I;
    end
    if (det(V) ~= 1)
        disp('U is not a rotations matrix, changing it to one by multiplication with matrix')
        V = V*I;
    end
    %% Berechnung der beiden möglichen Lösungen
    Rz1 = [0,-1,0;1,0,0;0,0,1];
    Rz2 = [0, 1, 0; -1, 0, 0; 0, 0, 1];
    
    % 2 Lösungne für R
    R1 = U*Rz1'*V';
    R2 = U*Rz2'*V';
    
    % 2 Lösungen für T
    T1 = U*Rz1*S*U';
    T2 = U*Rz2*S*U';
    
    T1 = [T1(3,2), T1(1,3), T1(2,1)]';
    T2 = [T2(3,2), T2(1,3), T2(2,1)]';

    
end
