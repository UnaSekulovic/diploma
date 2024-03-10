function [U, S, V] = tSVD(M)
    % Input: M ∈ R^(n_1 x n_2...x n_N)
    % rho = n_3*n_4*...*n_N
    
    % Extract dimensions
    sizes = size(M);
    N = numel(sizes);
    rho = prod(sizes(3:end));
    
    % Step 1: Apply FFT along dimensions 3 to N
    for i = 3:N
        D = fft(M, [], i);
    end
    
   
    U_1 = zeros(sizes);
    S_1 = zeros(sizes);
    V_1 = zeros(sizes);
    
    for i = 1:rho
        
        [U, S, V] = svd(D(:, :, i), 'econ');
        
        % Store results
        U_1(:, :, i) = U;
        S_1(:, :, i) = S;
        V_1(:, :, i) = V;
    end
    
    % Step 3: Apply IFFT along dimensions 3 to ρ
    for i = 3:rho
        U = ifft(U_1, [], i);
        S = ifft(S_1, [], i);
        V = ifft(V_1, [], i);
    end
    
    [U, S, V];
end
