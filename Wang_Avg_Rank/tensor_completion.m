function X = tensor_completion(M, omega, opts)
   
    [n1, n2, n3] = size(M);
    %X = P_Omega(M, omega);
    X = zeros(n1, n2, n3);
    Y = zeros(n1, n2, n3);
    
    % Set default options
    K = 200;
    rho = n3;
    mi = 0.5;
    lambda_min = 1;
    lambda_max = 5;
    delta = 1;
    epsilon = 1;
    
 
    if nargin > 2 && ~isempty(opts)
        if isfield(opts, 'max_iter');    K = opts.max_iter;    end
        if isfield(opts, 'rho');         rho = opts.rho;              end
    end
    
    lambda = lambda_max;
    
    for k = 1:K
        X_old = X;
        Y_old = Y;
        lambda_old = lambda;
        
        X_temp = (1 - mi)*Y_old + mi*X_old;
        threshold = lambda_old * mi;
        Y = THT(X_temp, threshold, n3);
        
        X = M .* omega + Y .*(1 - omega);
        
        norm_Y = norm(Y(:) - Y_old(:), Inf);
        norm_X = norm(X(:) - X_old(:), Inf);
        
        if min(norm_Y, norm_X) <= delta
            lambda = max(rho * lambda_old, lambda_min);
        end
        
        if min(norm_Y, norm_X) <= epsilon && min(norm_Y, norm_X) ~= 0
            break;
        end
    end
end


function D = THT(X, threshold, n3)
    for i = 1:n3
        [U, S, V] = svd(X(:, :, i), 'econ');
        S_threshold = max(S - threshold, 0);
        D(:, :, i) = U * S_threshold * V';
    end
end




