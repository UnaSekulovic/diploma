function rse_value = rse(X, X_hat)
    norm_X = norm(X, 'fro');
    %norm_X
    rse_value = norm(X(:) - X_hat(:), 'fro') / norm_X;
    %rse_value
end
