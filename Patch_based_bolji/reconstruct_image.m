
image = double(imread('giant.png'));
image = image / 255;
maxP = max(abs(image(:)));

patch_size = 32;
stride = 16;

figure(1);
imshow(image);
title('Originalna slika');

% zašumljivanje 
p = 0.8;
[n1, n2, n3] = size(image);
omega = find(rand(n1 * n2 * n3, 1) < p);
omega2 = zeros(n1, n2, n3);
Iones = ones(n1, n2, n3);
omega2(omega) = Iones(omega);

noisy_image = image .* omega2;

figure(2);
imshow(noisy_image);
title('Zašumljena slika');

% razdeliti sliko na patcheve
[patches, positions] = extract_patch(image, patch_size, stride);

% transformacija patcheva v vektorje, patch je bil velikosti 32x32x3, in
% zdaj ga transformiram v velikost (32^2)x3, ne vem ali je to prava
% velikost? ali sem prav razumela?
patch_vectors = cellfun(@(patch) patch_to_vector(patch), patches, 'UniformOutput', false);

% naključna izbira ciljnih patch-eva
num_target_patches = 10;
num_similar_patches = 30;
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
for i = 1:num_target_patches
    % rekonstrukcija vsake grupe
    group = similar_patches{i};
    %reconstructed_group = reconstruct_group(group, omega2);
    %patch_indices = sorted_indices(1:num_similar_patches);
    positions_for_group = similar_positions{i}; % Pozicije za trenutnu grupu
    reconstructed_group = reconstruct_group(group, omega2, positions_for_group, patch_size);
    
    % posodobiti rekonstruirane grupe
    for j = 1:num_similar_patches
        patch = reshape(reconstructed_group(:, j, :), [patch_size, patch_size, n3]);
        
        % iskanje pozicij patcheva
        % [row_idx, col_idx] = ind2sub(size(image) - patch_size + 1, target_indices(i));
        % %[row_idx, col_idx] = ind2sub(size(image) - patch_size + 1, patch_indices(j));
        % row_range = row_idx:row_idx + patch_size - 1;
        % col_range = col_idx:col_idx + patch_size - 1;

        pos = positions_for_group(j, :);
        row_range = pos(1):pos(1) + patch_size - 1;
        col_range = pos(2):pos(2) + patch_size - 1;
        
        reconstructed_patches(row_range, col_range, :) = reconstructed_patches(row_range, col_range, :) + patch;
        weight_matrix(row_range, col_range, :) = weight_matrix(row_range, col_range, :) + 1;
    end
end

reconstructed_image = reconstructed_patches ./ weight_matrix;

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
    
    opts.DEBUG = 1;
    reconstructed_group_tensor = Wang(group .* omega, omega, opts);
    
    reconstructed_group = reconstructed_group_tensor;
end
