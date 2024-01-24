function X = HaLRTC(T,omega)
    % Input: X_omega = T_omega, ro, K 
    % Output: X - Completed Tensor

    sz = size(T);
    [n1,n2,n3] = size(T);

    X = zeros(sz);
    X(omega) = T(omega);
    %X(omega) = M(omega);
    Y = zeros(sz);
    %X_nomega = 0
    
    K = 500;
    rho = 1.1;
    tau = 5*((n1+n2)/2);
    

    for k = 1:K
        tabela_M = zeros(numel(sz));
        % Iterate over each mode
        for i = 1:numel(sz)
            % Mode folding
            M = fold(D_tau(X(:, :, i) + (1/rho) * Y(:, :, i), tau), sz, i);
            tabela_M(i) = M;
        end
        
        X(omega) = (1/numel(sz)) * (sum(tabela_M) - (1/rho) * Y(i));

       tabela_Y = zeros(len(tabela_M));
            % Unfold the tensor along the current mode and update Y_i
        for i = 1:len(tabela_M)
            Y = Y - rho * (tabela_M(i) - X);
            tabela_Y(i) = Y;
        end
    end

    X;
end



function D_tau_result = D_tau(X, tau)
    
    [U, S, V] = svd(X);

    % Modify singular values according to the parameter tau
    S_tau = soft_thresholding(S, tau);

    % Reconstruct the tensor using modified singular values
    D_tau_result = U * S_tau * V';

end

function S_tau = soft_thresholding(S, tau)
    % Function for soft thresholding of singular values
    % Input: S - Singular values, tau - threshold parameter
    % Output: S_tau - Modified singular values

    S_tau = sign(S) .* max(0, abs(S) - tau);
    %sign ovdje ne treba!!
end

% function foldedTensor = foldi(T, mode)
%     % Function to perform mode folding along the specified mode
%     % Input: T - Tensor, mode - Mode along which folding is done
%     % Output: foldedTensor - Folded tensor
% 
%     sz = size(T);
% 
%     % Rearrange dimensions to fold along the specified mode
%     order = [mode, 1:mode-1, mode+1:length(sz)];
%     foldedTensor = permute(T, order);
%     foldedTensor = reshape(foldedTensor, sz(mode), []);
% 
%     % Transpose to match the order of dimensions in the pseudocode
%     foldedTensor = foldedTensor';
% 
% end
