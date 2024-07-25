
image = double(imread('0020.png'));
image = image / 255;
maxP = max(abs(image(:)));

patch_size = 32;
stride = 16;

%figure(1);
%imshow(image);
%title('Originalna slika');

% zašumljivanje 
p = 0.8;
[n1, n2, n3] = size(image);
omega = find(rand(n1 * n2 * n3, 1) < p);
omega2 = zeros(n1, n2, n3);
Iones = ones(n1, n2, n3);
omega2(omega) = Iones(omega);

noisy_image = image .* omega2;

%figure(2);
%imshow(noisy_image);
%title('Zašumljena slika');

% razdeliti sliko na patcheve
[patches, positions] = extract_patch(image, patch_size, stride);

% transformacija patcheva v vektorje, patch je bil velikosti 32x32x3, in
% zdaj ga transformiram v velikost (32^2)x3, ne vem ali je to prava
% velikost? ali sem prav razumela?
patch_vectors = cellfun(@(patch) patch_to_vector(patch), patches, 'UniformOutput', false);

% naključna izbira ciljnih patch-eva
num_target_patches = 90;
num_similar_patches = 10;
target_indices = randperm(length(patch_vectors), num_target_patches);
target_patches = patch_vectors(target_indices);


% iskanje podobnih patchev za vsak ciljni, tako da je potem velikost vsake
% grupe (tenzorja) (32^2)x30x3
similar_patches = cell(num_target_patches, 1);
similar_positions = cell(num_target_patches, 1);
for i = 1:num_target_patches
    target_patch = target_patches{i};
    distances = cellfun(@(patch) sum((patch - target_patch).^2, 'all'), patch_vectors);
    [~, sorted_indices] = sort(distances);
    
    group = zeros(patch_size^2, num_similar_patches, n3);
    
    for j = 1:num_similar_patches
        similar_patch_vector = patch_vectors{sorted_indices(j)};
        for c = 1:n3
            group(:, j, c) = similar_patch_vector(:, c);
        end
    end
    
    similar_patches{i} = group;
    similar_positions{i} = positions(sorted_indices(1:num_similar_patches), :);
end

% rekonstrukcija
reconstructed_patches = zeros(size(image));
weight_matrix = zeros(size(image));
tic;
for i = 1:num_target_patches
    % rekonstrukcija vsake grupe
    group = similar_patches{i};
    positions_for_group = similar_positions{i}; 
    reconstructed_group = reconstruct_group(group, omega2, positions_for_group, patch_size);
    
    % posodobiti rekonstruirane grupe
    for j = 1:num_similar_patches
        patch = reshape(reconstructed_group(:, j, :), [patch_size, patch_size, n3]);

        pos = positions_for_group(j, :);
        row_range = pos(1):pos(1) + patch_size - 1;
        col_range = pos(2):pos(2) + patch_size - 1;
        
        reconstructed_patches(row_range, col_range, :) = reconstructed_patches(row_range, col_range, :) + patch;
        weight_matrix(row_range, col_range, :) = weight_matrix(row_range, col_range, :) + 1;
    end
end

%weight_matrix(weight_matrix == 0) = 1;
reconstructed_image = reconstructed_patches ./ weight_matrix;

unreconstructed_mask = (weight_matrix == 0);
reconstructed_image(unreconstructed_mask) = noisy_image(unreconstructed_mask);

%TO je TA ZADNJI DEL, KI JE ZDAJ DODAN
% iskanje patchev, ki niso bili rekonstruirani
unreconstructed_patches = extract_unreconstructed_patches(reconstructed_image, unreconstructed_mask, patch_size, stride);

% rekonstrukcija nerekonstruiranih patchev
for k = 1:length(unreconstructed_patches)
    patch_info = unreconstructed_patches{k};
    patch_data = patch_info.data;
    patch_mask = patch_info.mask;

    opts.DEBUG = 1;
    reconstructed_patch = lrtc_tnn(patch_data .* patch_mask, patch_mask, opts);

    row_range = patch_info.row_range;
    col_range = patch_info.col_range;

    reconstructed_image(row_range, col_range, :) = reconstructed_patch;
end

kraj = toc;
psnr = PSNR(image,reconstructed_image,maxP);

rse1 = rse(image, reconstructed_image);

%saveas(imshow(reconstructed_image), '0020_tnn_90_10_druga.png');
fprintf('PSNR: %.4f, num_target: %.4f, num_group: %.4f\n', psnr, num_target_patches, num_similar_patches);
fprintf('RSE: %.4f, num_target: %.4f, num_group: %.4f\n', rse1, num_target_patches, num_similar_patches);
fprintf('Time: %.4f, num_target: %.4f, num_group: %.4f\n', kraj, num_target_patches, num_similar_patches);
                


figure(3);
subplot(1, 2, 1);
imshow(noisy_image);
title('Zašumljena slika');

subplot(1, 2, 2);
imshow(reconstructed_image);
%imshow(uint8(reconstructed_image));
title('Rekonstruisana slika');


function [patches, positions] = extract_patch(image, patch_size, stride)
    [n1, n2, ~] = size(image);
    patches = {};
    positions = [];
    st = 1;
    for i = 1:stride:n1 - patch_size + 1
        for j = 1:stride:n2 - patch_size + 1
            patch = image(i : i + patch_size - 1, j : j + patch_size - 1, :);
            patches{st} = patch;
            positions(st, :) = [i, j];
            st = st + 1;
        end
    end
end

function patch_vector = patch_to_vector(patch)
    [patch_size, ~, n3] = size(patch);
    
    patch_vector = zeros(patch_size^2, n3);
    
    for c = 1:n3
        patch_channel = patch(:, :, c);
        patch_vector(:, c) = patch_channel(:);
    end
end


function reconstructed_group = reconstruct_group(group, global_omega, positions_for_group, patch_size)
    [~, num_patches, n3] = size(group);
    
    omega = zeros(size(group));
    for idx = 1:num_patches
        %[row, col] = ind2sub([size(global_omega, 1) - patch_size + 1, size(global_omega, 2) - patch_size + 1], patch_indices(idx));
        pos = positions_for_group(idx, :);
        for c = 1:n3
            %patch_omega = global_omega(row:row+patch_size-1, col:col+patch_size-1, c);
            patch_omega = global_omega(pos(1):pos(1)+patch_size-1, pos(2):pos(2)+patch_size-1, c);
            omega(:, idx, c) = patch_omega(:);
        end
    end
    omega2 = omega;
    omega2 = omega2 == 1;
    opts.DEBUG = 1;
    reconstructed_group_tensor = lrtc_tnn(group .* omega2, omega2, opts);
    
    reconstructed_group = reconstructed_group_tensor;
end

function unreconstructed_patches = extract_unreconstructed_patches(image, mask, patch_size, stride)
    [n1, n2, ~] = size(image);
    unreconstructed_patches = {};
    
    st = 1;
    for i = 1:stride:n1 - patch_size + 1
        for j = 1:stride:n2 - patch_size + 1
            patch_mask = mask(i:i+patch_size-1, j:j+patch_size-1, :);
            if any(patch_mask(:))
                patch_data = image(i:i+patch_size-1, j:j+patch_size-1, :);
                unreconstructed_patches{st} = struct('data', patch_data, 'mask', patch_mask, 'row_range', i:i+patch_size-1, 'col_range', j:j+patch_size-1);
                st = st + 1;
            end
        end
    end
end