function X_result = tensor_completion(M, Omega, opts)
    [n1, n2, n3] = size(M);
    X = P_Omega(M, Omega);
    Z = zeros(n1, n2, n3);
    Q = zeros(n1, n2, n3);
    sizes = size(M);
    
    % Set default options
    
    %vrednosti za 'tol' sem vzela iz drugega algoritma (TNN), ki je že
    %implementiran
    tol = 1e-6; 
    %to vrednost sem dala naključno
    max_iter = 200;
    %to sem dala kot v algoritmu tSVD na 3. strani (na samem začetku
    %strane)
    rho = prod(sizes(3:end));
    %to tudi ne vem na koju vrednost naj nastavim
    threshold = 5*((n1+n2)/2); 
    
    %to za opts sem tudi samo vzela iz tega implementiranega algoritma
    % Update options if provided
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
        
        % Update X
        %Y = P_Omega(X_old - Z_old + Q_old, Omega);
        
        %
        %Y = P_Omega(X_old, Omega);
        
        %formula 8 na strani 5 in stran 3, odstavek 2
        X = update_X(X_old, Z_old, Q_old);
        
        %strana 5 - formula 9, ta del s frobeniusovo normo sem naredila
        %enako kot za 'update_X' (ne vem ce je to prav), del s
        %nuklearno normo - strana 3, Theorem 2.4.1.
        Z = update_Z(X, Z_old, Q_old, rho, threshold);
        
        % strana 5, formula 10
        Q = Q_old + (X - Z);
        
        % tudi to preverjanje pogoja sem samo vzela iz že implemetiranega
        % algoritma, ker ne vem kaj je pogoj tukaj
        if max(norm(X(:) - X_old(:)), norm(Z(:) - Z_old(:))) < tol
            break;
        end
    end
    
    X_result = X;
end


function X = update_X(Y, Z, Q)
   
    [U, S, V] = tSVD(Y - Z + Q);

    %ne vem kako po tSVD dobiti X!!!
    tempV(:,:,1)=V(:,:,1)';
    tempV(:,:,2)=V(:,:,3)';
    tempV(:,:,3)=V(:,:,2)';
    X = tprod(tprod(U,S),tempV);

   
    %X satisfies constraint P_omega(X) = P_omega(M)
    %ne vem ce se to da izraziti tko!!!
    %size(X)
    %size(Omega)
    %pause
    %X
    %pause
    %X = X .* Omega + M .* (1 - Omega);
    %X
    %pause   
end


function Z = update_Z(X, Z, Q, rho, threshold)
    % Transform to Fourier domain along the third mode - ker je napisano X
    % stresica, in enako za Z i Q - strana 5, odstavek 3
     X_fft = fft(X, [], 3);
     Z_fft = fft(Z, [], 3);
     Q_fft = fft(Q, [], 3);
     
     
    %prvi del s nuklearno normo
    
    [U, S, V] = tSVD(Z_fft);
    
    
    S_threshold = max(S - threshold, 0);
    
    tempV(:,:,1)=V(:,:,1)';
    tempV(:,:,2)=V(:,:,3)';
    tempV(:,:,3)=V(:,:,2)';
    Z_new = tprod(tprod(U,S_threshold),tempV);
    
    nuclear = sum(Z_new);
    
    %drugi del - frobenius 
    
    [U, S, V] = tSVD(Z_fft - (X_fft + Q_fft));
    tempV(:,:,1)=V(:,:,1)';
    tempV(:,:,2)=V(:,:,3)';
    tempV(:,:,3)=V(:,:,2)';
    Z_fro = tprod(tprod(U,S),tempV);
    
    Z_temp = 1/rho * nuclear + 0.5 * Z_fro;

    % Perform the inverse Fourier transform to get Z
    Z = ifftn(Z_temp);

  
end


%projection of known values
function Y = P_Omega(X, omega)
    Y = X .* omega;
end

% function Y = tTranspose(X)
%     % Transpose tensor X along the third mode
%     Y = permute(X, [1, 3, 2]);
% end













