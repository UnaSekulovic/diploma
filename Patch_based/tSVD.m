function X_result = tensor_completion(M, Omega, opts)
    [n1, n2, n3] = size(M);
    X = zeros(n1, n2, n3);
    Z = zeros(n1, n2, n3);
    Q = zeros(n1, n2, n3);
    
    
    % Set default options
    tol = 1e-6; 
    max_iter = 200;
    rho = n3;
 
    if nargin > 2 && ~isempty(opts)
        if isfield(opts, 'tol');         tol = opts.tol;              end
        if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
        if isfield(opts, 'rho');         rho = opts.rho;              end
    end
    
    % Main loop
    for iter = 1:max_iter
        X_old = X;
        Z_old = Z;
        Q_old = Q;
        
        %formula 8 na strani 5 in stran 3, odstavek 2
        %X = update_X(X_old, Z_old, Q_old, M, Omega);
        X = M .* Omega + (Z_old-Q_old) .* (1 - Omega);
        
        %strana 5 - formula 9, ta del s frobeniusovo normo sem naredila
        %enako kot za 'update_X' (ne vem ce je to prav), del s
        %nuklearno normo - strana 3, Theorem 2.4.1.
        Z = update_Z(X, Z_old, Q_old, rho);
        %norm(Z(:)-Z_old(:))
        
        % strana 5, formula 10
        Q = Q_old + (X - Z);
        
        % tudi to preverjanje pogoja sem samo vzela iz že implemetiranega
        % algoritma, ker ne vem kaj je pogoj tukaj
       % norm(X(:) - X_old(:))
        if max(norm(X(:) - X_old(:)), norm(Z(:) - Z_old(:))) < tol
            break;
        end
    end
    
    X_result = X;
end





function Z = update_Z(X, Z, Q, rho)
    % Transform to Fourier domain along the third mode - ker je napisano X
    % stresica, in enako za Z i Q - strana 5, odstavek 3
     X_fft = fft(X, [], 3);
     Q_fft = fft(Q, [], 3);
    
    
    [U1, S1, V1] = svd(X_fft(:, :, 1) + Q_fft(:, :, 1));
    [U2, S2, V2] = svd(X_fft(:, :, 2) + Q_fft(:, :, 2));
    [U3, S3, V3] = svd(X_fft(:, :, 3) + Q_fft(:, :, 3));
    %S(1, 1)
    
    threshold = 1 / rho;
    S_threshold1 = max(S1 - threshold, 0);
    S_threshold2 = max(S2 - threshold, 0);
    S_threshold3 = max(S3 - threshold, 0);
    

    Z(:, :, 1) = U1 * S_threshold1 * V1';
    Z(:, :, 2) = U2 * S_threshold2 * V2';
    Z(:, :, 3) = U3 * S_threshold3 * V3';
    
    Z = ifft(Z, [], 3);

  
end















