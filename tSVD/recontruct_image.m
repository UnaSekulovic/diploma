%clear
M = double(imread('giant.png'));
%M = M;

%M

X = tSVD(M);

%X

figure(1)
%imshow(X, [])

imagesc(X);
%colormap('hot');
%colorbar;  % Display colorbar for reference
caxis([min(X(:)), max(X(:))]);

% Example usage:
% Assuming you have a tensor M that you want to process with tSVD
% For demonstration purposes, you can use a random tensor M:
%M = randn(100, 100, 5);  % Replace this with your actual data

% Apply tSVD
%reconstructed_image = tSVD(M);

% Visualize the reconstructed image
%imshow(reconstructed_image, []);

% You may also consider using imagesc instead of imshow, depending on your data
% imagesc(reconstructed_image);
