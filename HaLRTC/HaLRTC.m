function X = HaLRTC(T, omega, opts)
   %M_i - razumela sem da obstaja več tenzorjev, kot M_1, M_2 itn
   %enako za Y_i
   %ker 'i' gre od 1 do n, jaz sem razumela da je 'n' enako številu dimenzij
   %tenzora (n-mode tensor), in v našem primeru je to enako 3, tko da sem jaz takoj na
   %začetku ročno dala M1, M2, M3, in enako za Y
   
   %X_(i) - ko je indeks napisen s oklepaji to sem razumela da pomeni unfold
   %po tej dimenziji, to sem tako razumela zato ker indeks s oklepaji sem
   %najdla na strani 2, poglavje 1.1 Notation, desni stolpec, 2. odstavek,
   %kjer je razloženo kaj je unfold = X_(k)


    [n1, n2, n3] = size(T);
    X = P_Omega(T, omega);
    Y1 = zeros(n1, n2, n3);
    Y2 = zeros(n1, n2, n3);
    Y3 = zeros(n1, n2, n3);
    M1 = zeros(n1, n2, n3);
    M2 = zeros(n1, n2, n3);
    M3 = zeros(n1, n2, n3);
    
    
    
    % Set default options
    %tol = 1e-6; 
    K = 10;
    rho = n3;
    threshold = 1 / rho;
 
    if nargin > 2 && ~isempty(opts)
        %if isfield(opts, 'tol');         tol = opts.tol;              end
        if isfield(opts, 'max_iter');    K = opts.max_iter;    end
        if isfield(opts, 'rho');         rho = opts.rho;              end
    end
    
    
    for k = 1:K
        for i = 1:3
            if i == 1
                %"LALA"
                temp1 = unfold(X, i) + 1/rho * unfold(Y1, i);
                [U, S, V] = svd(temp1);
                S_threshold = max(S - threshold, 0);
                M1 = fold(U * S_threshold * V', i, [n1, n2, n3]);
                %M1
            elseif i == 2
                %"LA2"
                temp2 = unfold(X, i) + 1/rho * unfold(Y2, i);
                [U, S, V] = svd(temp2);
                S_threshold = max(S - threshold, 0);
                M2 = fold(U * S_threshold * V', i, [n1, n2, n3]);
               % M2
            elseif i == 3
                %"LA3"
                temp3 = unfold(X, i) + 1/rho * unfold(Y3, i);
                [U, S, V] = svd(temp3);
                S_threshold = max(S - threshold, 0);
                M3 = fold(U * S_threshold * V', i, [n1, n2, n3]);
            end
        end
        
        X = 1/3 * ((M1 + M2 + M3) - 1/rho * Y1);
        
        Y1 = Y1 - rho * (M1 - X);
        Y2 = Y2 - rho * (M2 - X);
        Y3 = Y3 - rho * (M3 - X);
    end
end

%projection of known values
function Y = P_Omega(X, omega)
    Y = X .* omega;
end




