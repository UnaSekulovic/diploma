%prebrati sliko
clear
image = double(imread('giant.png'));
image = image/255; 
maxP = max(abs(image(:)));
%size(image)
patch_size = 250; 
stride = 4; 

figure(1);
imshow(image);
%zašumiti sliko
p = 0.5; 
[n1,n2,n3] = size(image);
omega = find(rand(n1*n2*n3,1)<p);
omega2 = zeros(n1,n2,n3);
Iones = ones(n1,n2,n3);
omega2(omega) = Iones(omega);


%razdeliti sliko na patcheve
patches = extract_patch(image, patch_size, stride);
%size(patches, 2)

%rekonstrukcija patcheva
reconstructed = {};
for i = 1 : size(patches, 2)
    reconstructed{i} = reconstruct_patch(patches{i}, omega2);
end

reconstructed_image = zeros(size(image));
count = zeros(size(image));
br = 1;
for i = 1 : stride : size(image, 1) - patch_size + 1
    for j = 1 : stride : size(image, 2) - patch_size + 1
        reconstructed_image(i : i + patch_size - 1, j : j + patch_size - 1, :) = reconstructed_image(i : i + patch_size - 1, j : j + patch_size - 1, :) + reconstructed{br};
        count(i : i + patch_size - 1, j : j + patch_size - 1, :) = count(i : i + patch_size - 1, j : j + patch_size - 1, :) + 1;
        br = br + 1;
    end
end

%br
reconstructed_image = reconstructed_image ./ count;

figure(2)
%imshow(uint8(reconstructed_image));
imshow(reconstructed_image);

psnr = PSNR(image, reconstructed_image, maxP)

function patches = extract_patch(X, patch_size, stride)
    [n1, n2, ~] = size(X);
    patches = {};
    st = 1;
    for i = 1 : stride : n1 - patch_size + 1
        for j = 1 : stride : n2 - patch_size + 1
            patch = X(i : i + patch_size - 1, j : j + patch_size - 1, :);
            patches{st} = patch;
            st = st + 1;
        end
    end
end



%funkcija za rekonstrukcijo patcheva
function reconstructed = reconstruct_patch(patch, global_omega)
    [n1,n2,n3] = size(patch);

    omega = zeros(n1, n2, n3);
    for i = 1:n1
        for j = 1:n2
            for k = 1:n3
                if global_omega(i, j, k) == 1
                    omega(i, j, k) = global_omega(i, j, k);
                end
            end
        end
    end

    M = zeros(n1,n2,n3);
    M(omega == 1) = patch(omega == 1);

    omega2 = zeros(n1,n2,n3);
    omega2(omega == 1) = 1;

    opts.DEBUG = 1;
    
    reconstructed = Wang(M, omega2, opts);
end

