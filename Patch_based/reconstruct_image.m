% clear
% X = double(imread('giant.png'));
% X = X/255;    
% maxP = max(abs(X(:)));
% [n1,n2,n3] = size(X);
% 
% 
% figure(1)
% imshow(X)
% 
% p = 0.5; % sampling rate
% 
% omega = find(rand(n1*n2*n3,1)<p);
% M = zeros(n1,n2,n3);
% M(omega) = X(omega);
% 
% omega2 = zeros(n1,n2,n3);
% Iones = ones(n1,n2,n3);
% omega2(omega) = Iones(omega);
% 
% figure(2)
% imshow(M)
% opts.DEBUG = 1;
% 
% X_rez = tensor_completion(M,omega2,opts);
% psnr = PSNR(X,X_rez,maxP)
% 
% 
% figure(3)
% imshow(X_rez)



clear
image = imread('giant.png');
patch_size = 8; 
stride = 4; 

%razdeliti sliko na patcheve
patches = extract_patch(image, patch_size, stride);

%rekonstrukcija patcheva
reconstructed = zeros(size(patches));
for i = 1 : size(patches, 4)
    reconstructed(:, :, :, i) = reconstruct_patch(patches(:, :, :, i));
end

%spajanje patcheva
reconstructed_image = zeros(size(image));
count = zeros(size(image));
br = 1;
for i = 1 : stride : size(image, 1) - patch_size + 1
    for j = 1 : stride : size(image, 2) - patch_size + 1
        reconstructed_image(i : i + patch_size - 1, j : j + patch_size - 1, :) = reconstructed_image(i : i + patch_size - 1, j : j + patch_size - 1, :) + reconstructed(:, :, :, br);
        count(i : i + patch_size - 1, j : j + patch_size - 1, :) = count(i : i + patch_size - 1, j : j + patch_size - 1, :) + 1;
        br = br + 1;
    end
end
reconstructed_image = reconstructed_image ./ count;


imshow(uint8(reconstructed_image));


function patches = extract_patch(X, patch_size, stride)
    [n1, n2, ~] = size(X);
    patches = [];
    for i = 1 : stride : n1 - patch_size + 1
        for j = 1 : stride : n2 - patch_size + 1
            patch = X(i : i + patch_size - 1, j : j + patch_size - 1, :);
           % patch
            patches = cat(4, patches, patch);
        end
    end
end



%funkcija za rekonstrukcijo patcheva
function reconstructed = reconstruct_patch(patch)
    p = 0.5; 
    [n1,n2,n3] = size(patch);
    omega = find(rand(n1*n2*n3,1)<p);
    M = zeros(n1,n2,n3);
    M(omega) = patch(omega);
 
    omega2 = zeros(n1,n2,n3);
    Iones = ones(n1,n2,n3);
    omega2(omega) = Iones(omega);

    opts.DEBUG = 1;
    
    reconstructed = tensor_completion(M, omega2, opts);
end

