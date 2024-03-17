% %clear
% M = double(imread('giant_missing_60.png'));
% %M = imread('lena_missing_60.png');
% %M 
% 
% %M
% sizes = size(M);
% N = numel(sizes);
% rho = prod(sizes(3:end));
% %X = tensor_completion(M, M, 100, rho);
% [U, S, V] = tSVD(M);
% 
% %X
% 
% figure(1)
% %imshow(X, [])
% 
% imagesc(X);



% applying tnsor completion for image inpainting

clear
X = double(imread('giant.png'));
X = X/255;    
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);

p = 0.5; % sampling rate

omega = find(rand(n1*n2*n3,1)<p);
M = zeros(n1,n2,n3);
M(omega) = X(omega);
%M
%[n1, n2, n3] = size(M)

%M2 = Frontal2Lateral(M); % each lateral slice is a channel of the image
omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
omega2(omega) = Iones(omega);
%omega2 = Frontal2Lateral(omega2);
%omega2 = find(omega2==1);

opts.DEBUG = 1;
%omega2
X = tensor_completion(M,omega2,opts);
%Xhat = tensor_completion(M, omega2);
%Xhat = max(Xhat,0);
%Xhat = min(Xhat,maxP);
%Xhat = Lateral2Frontal(Xhat); % each lateral slice is a channel of the image
%psnr = PSNR(X,Xhat,maxP)


figure(1)
% subplot(1,3,1)
% imshow(X)
% subplot(1,3,2)
% imshow(M)
% subplot(1,3,3)
imshow(X)


