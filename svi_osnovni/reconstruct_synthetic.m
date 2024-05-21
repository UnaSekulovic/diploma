clear

%imena_algoritama = {'HaLRTC', 'tSVD', 'Wang'};
%fid = fopen('rezultati_syntetic_TNN.txt', 'w');

%for j = 1:length(imena_algoritama)
    %algoritam = str2func(imena_algoritama{j});
    %for r=1:1:30
    r = 2;
        %for p = 0.01 : 0.05: 0.99
        p = 0.8;
            psnr_all = [];
            time = [];
            rse_all = [];
            for i = 1 : 10
                P = rand(481, r, 3);
                Q = rand(r, 321, 3);
                %X = P .* Q;
                X = tprod(P, Q);
   
                maxP = max(abs(X(:)));
                [n1,n2,n3] = size(X);
                omega = find(rand(n1*n2*n3,1)<p);
                M = zeros(n1,n2,n3);
                M(omega) = X(omega);

                omega2 = zeros(n1,n2,n3);
                Iones = ones(n1,n2,n3);
                omega2(omega) = Iones(omega);
                omega3 = find(omega2==1);
                opts.DEBUG = 1;

                tic;
                %if strcmp(imena_algoritama{j}, 'lrtc_tnn')
                    %X_rez = lrtc_tnn(M,omega3,opts);
                    X_rez = tSVD(M, omega2, opts);
                %else
                 %   X_rez = algoritam(M,omega2,opts);
                %end

                kraj = toc;
        
                psnr = PSNR(X,X_rez,maxP);
                
                psnr_all(i) = psnr;
                time(i) = kraj;
                rse_all(i) = rse(X, X_rez);
            end
            %povprečni PSNR
            psnr_avg = mean(psnr_all);
            fprintf('Average PSNR - tSVD: %.4f za p: %d in r: %d\n', psnr_avg, p, r);
            %fprintf(fid, 'Average PSNR - TNN: %.2f za p: %d in r: %d\n', psnr_avg, p, r);

            %čas
            cas = mean(time);
            fprintf('Time - tSVD: %.4f za p: %d in r: %d\n', cas, p, r);
            %fprintf(fid, 'Time - TNN: %.2f za p: %d in r: %d\n', cas, p, r);

            %povprečni RSE
            rse_avg = mean(rse_all);
            fprintf('Average RSE - tSVD: %.4f za p: %d in r: %d\n', rse_avg, p, r);
            %fprintf(fid, 'Average RSE - TNN: %.2f za p: %d in r: %d\n', rse_avg, p, r);
        %end
    %end
    
%end

%fclose(fid);


