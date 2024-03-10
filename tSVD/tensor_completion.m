function X_result = tensor_completion(M, Omega, opts)
    % Initialize variables
    [n1, n2, n3] = size(M);
    X = zeros(n1, n2, n3);
    Z = zeros(n1, n2, n3);
    Q = zeros(n1, n2, n3);
    %iter = 0;
    
    % Set default options
    tol = 1e-6; 
    max_iter = 500;
    %rho = 1.1;
    %mu = 1e-4;
    %max_mu = 1e10;
    %DEBUG = 0;
    
    % Update options if provided
    if nargin > 2 && ~isempty(opts)
        if isfield(opts, 'tol');         tol = opts.tol;              end
        if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
        %if isfield(opts, 'rho');         rho = opts.rho;              end
        %if isfield(opts, 'mu');          mu = opts.mu;                end
        %if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
        %if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end
    end
    
    % Main loop
    for iter = 1:max_iter
        X_old = X;
        Z_old = Z;
        
        % Update X
        Y = P_Omega(X - Z + Q, Omega);
        X = update_X(Y, M, Z, Q, Omega);
        
        % Update Z
        Z = update_Z(X, Z, Q, Omega, M);
        
        % Update Lagrange multiplier Q
        Q = Q + (X - Z);
        
        % Check convergence
        if max(norm(X(:) - X_old(:)), norm(Z(:) - Z_old(:))) < tol
            break;
        end
    end
    
    X_result = X;
end


function X = update_X(Y, M, Z, Q, Omega)
    % Transform to Fourier domain along the third mode
    %Y_hat = fft(Y, [], 3);
    %M_hat = fft(M, [], 3);
    %Z_hat = fft(Z, [], 3);
    %Q_hat = fft(Q, [], 3);

    %tSVD
    [U, S, V] = tSVD(Y - Z + Q);

%     % Transform back to spatial domain
%     U = ifft(U_hat, [], 3);
%     S = ifft(S_hat, [], 3);
%     V = ifft(V_hat, [], 3);
%     
%     %[n1, n2, n3] = size(U)
%     %[n1, n2, n3] = size(S)
%     %[n1, n2, n3] = size(V)
% 
%     % Update X based on the modified factors
%     X = tprod(tprod(U, S), tTranspose(V));
%     %X = U .* S .* V;
%     %X = tprod(U, [1 3], S, [1 2], V, [3 2]);

    %ne vem kako po tSVD dobiti X!!!
    X = U .* S .* conj(V);

    % Transform back to spatial domain
    %X = ifft(X_hat, [], 3);
    
    
     %X = your_tSVD_algorithm(Y - Z + Q);

    %X satisfies constraint P_omega(X) = P_omega(M)
    %ne vem ce se to da izraziti tko!!!
    X = X .* Omega + M .* (1 - Omega);
end


function Z = update_Z(X, Z, Q, Omega, M)
    % Transform to Fourier domain along the third mode
%     X_hat = fft(X, [], 3);
%     Z_hat = fft(Z, [], 3);
%     Q_hat = fft(Q, [], 3);

    % Perform t-SVD on the transformed tensors
    [U, S, V] = tSVD(Z + Q - X);

    % Modify the obtained factors according to the constraints
    % (No modifications needed for t-SVD)

%     % Transform back to spatial domain
%     U = ifft(U_hat, [], 3);
%     S = ifft(S_hat, [], 3);
%     V = ifft(V_hat, [], 3);
% 
%     % Update Z based on the modified factors
%     Z = tprod(tprod(U, S), tTranspose(V));
%     %Z = U .* S .* V;

    %ne vem ce je to tko!!
    Z = U .* S .* conj(V);

    % Transform back to spatial domain
    %Z = ifft(Z_hat_new, [], 3);
    
    
    %Z = your_tSVD_algorithm(X + Q);

    % Perform singular value thresholding
    % (Modify this part according to your thresholding implementation)
    % Z = singular_value_thresholding(Z);

    %Z satisfies the constraint P_omega(X) = P_omega(M)
    Z = Z .* Omega + M .* (1 - Omega);
end


%projection of known values
function Y = P_Omega(X, omega)
    Y = X .* omega;
end

% function Y = tTranspose(X)
%     % Transpose tensor X along the third mode
%     Y = permute(X, [1, 3, 2]);
% end













