clear
X = double(imread('giant.png'));
X = X/255;    
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);


figure(1)
imshow(X)

p = 0.4; % sampling rate

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
omega2(omega) = Iones(omega);

figure(2)
imshow(M)
opts.DEBUG = 1;

X_rez = HaLRTC(M,omega2,opts);
psnr = PSNR(X,X_rez,maxP)


figure(3)
imshow(X_rez)


