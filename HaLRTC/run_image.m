clear
T = double(imread('lena_missing_60.png'));
T = T/255;    
%maxP = max(abs(X(:)));
[n1,n2,n3] = size(T);

p = 0.5; % sampling rate

omega = find(rand(n1*n2*n3,1)<p);
%M = zeros(n1,n2,n3);
%M(omega) = X(omega);

X = HaLRTC(T, omega);

figure(1)
%imshow(X, [])

imagesc(X);