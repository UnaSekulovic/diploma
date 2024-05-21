RSE = [0.95, 0.48, 0.19, 0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.98, 0.78, 0.64, 0.51, 0.39, 0.27, 0.15, 0.04, 0.00, 0.00, 0.98, 0.84, 0.73, 0.63, 0.53, 0.43, 0.33, 0.22, 0.10, 0.00, 0.99, 0.87, 0.78, 0.69, 0.60, 0.51, 0.41, 0.30, 0.19, 0.05, 0.99, 0.88, 0.80, 0.72, 0.64, 0.55, 0.46, 0.35, 0.24, 0.09, 0.99, 0.89, 0.82, 0.74, 0.66, 0.58, 0.49, 0.39, 0.27, 0.13, 0.99, 0.90, 0.83, 0.75, 0.68, 0.59, 0.51, 0.41, 0.30, 0.15, 0.99, 0.91, 0.83, 0.76, 0.69, 0.61, 0.52, 0.43, 0.32, 0.17, 0.99, 0.91, 0.84, 0.77, 0.70, 0.62, 0.54, 0.44, 0.33, 0.19, 0.99, 0.91, 0.85, 0.78, 0.70, 0.63, 0.54, 0.45, 0.34, 0.20];

threshold = 0.01;

br = 1;

p_values = 0.01:0.1:0.99;
r_values = 1:3:30;

color_map = zeros(length(r_values), length(p_values), 3);

i=1;
for r=1:3:30
    j=1;
    for p = 0.01 : 0.1: 0.99
        if RSE(br) < threshold
            color_map(i, j, :) = [1, 1, 0];
        else 
            color_map(i, j, :) = [0, 0, 1];
        end
        br = br + 1;
        j = j + 1;
    end
    i = i + 1;
end

figure;

imshow(color_map, 'InitialMagnification', 'fit');


axis on;


% Podesite oznake na osama (tick marks) i etikete (tick labels)
xticks(1:length(p_values));
xticklabels(arrayfun(@num2str, p_values, 'UniformOutput', false));
yticks(1:length(r_values));
yticklabels(arrayfun(@num2str, r_values, 'UniformOutput', false));

% Podešavanje oznaka i naslova grafika
xlabel('p');
ylabel('r');
title('Wang');

saveas(gcf, 'graphRSE_Wang.png');