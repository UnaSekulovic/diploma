% clear
% T = double(imread('lena_missing_60.png'));
% T = T/255;    
% %maxP = max(abs(X(:)));
% [n1,n2,n3] = size(T);
% 
% p = 0.5; % sampling rate
% 
% omega = find(rand(n1*n2*n3,1)<p);
% %M = zeros(n1,n2,n3);
% %M(omega) = X(omega);
% 
% X = HaLRTC(T, 50);
% 
% figure(1)
% %imshow(X, [])
% 
% imagesc(X);


clear
X = double(imread('giant.png'));
X = X/255;    
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);


% % Inicijalizacija tensora omega
% omega = zeros(size(X));
% 
% % Ručno označavanje prisutnih regiona slike
% imshow(X); % Prikazuje sliku radi ručnog označavanja
% h = drawfreehand; % Alat za ručno označavanje regiona slike
% roi = createMask(h); % Kreira masku označenih regiona
% omega(roi) = 1; % Postavlja vrednosti u tensoru omega na 1 za označene regione
% 
% % Konverzija u binarni tensor
% omega = logical(omega);


%figure(1)
%imshow(X)

p = 0.5; % sampling rate

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);

omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
omega2(omega) = Iones(omega);

%figure(2)
%imshow(M)
opts.DEBUG = 1;

X_rez = HaLRTC(M,omega2,opts);
psnr = PSNR(X,X_rez,maxP)


figure(3)
imshow(X)


