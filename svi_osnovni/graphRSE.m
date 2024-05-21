RSE = [0.95, 0.62, 0.41, 0.23, 0.09, 0.01, 0.00, 0.00, 0.00, 0.00, 0.98, 0.84, 0.74, 0.64, 0.54, 0.44, 0.34, 0.23, 0.11, 0.00, 0.99, 0.88, 0.80, 0.72, 0.63, 0.55, 0.45, 0.35, 0.24, 0.09, 0.99, 0.90, 0.82, 0.75, 0.67, 0.59, 0.50, 0.41, 0.29, 0.15, 0.99, 0.91, 0.84, 0.77, 0.69, 0.61, 0.53, 0.44, 0.33, 0.18, 0.99, 0.92, 0.85, 0.78, 0.71, 0.63, 0.55, 0.46, 0.35, 0.20, 0.99, 0.92, 0.85, 0.79, 0.72, 0.64, 0.56, 0.47, 0.36, 0.22, 0.99, 0.92, 0.86, 0.79, 0.72, 0.65, 0.58, 0.48, 0.37, 0.23, 0.99, 0.93, 0.86, 0.80, 0.73, 0.66, 0.57, 0.48, 0.38, 0.24, 0.99, 0.93, 0.86, 0.80, 0.73, 0.66, 0.58, 0.49, 0.38, 0.24];

threshold = 0.001;

br = 1;

p_values = 0.01:0.1:0.99;
r_values = 1:3:30;
% Inicijalizujte matricu za prikaz
%color_map = zeros(11, 10, 3);
color_map = zeros(length(r_values), length(p_values), 3);

i=1;
for r=1:3:30
    j=1;
    for p = 0.01 : 0.1: 0.99
        if RSE(br) < threshold
            %scatter(p, r, "yellow", 'filled'); % Uspešne rekonstrukcije su plave
            %scatter(p_values(j), r_values(i), 100, successful_color, 'filled');
            %successful_points = [successful_points; p, r];
            %color_map(i, j, :) = [0, 0, 1];
            color_map(i, j, :) = [1, 1, 0];
            %fprintf("%d \n", br);
        else 
            %scatter(p, r, "blue", 'filled'); % Uspešne rekonstrukcije su plave
            %unsuccessful_points = [unsuccessful_points; p, r];
            color_map(i, j, :) = [0, 0, 1];
            %fprintf("LALA");
            %fprintf("%d \n", RSE[br]);
        end
        br = br + 1;
        j = j + 1;
    end
    i = i + 1;
end

figure;
% Nacrtajte graf koristeći imshow
imshow(color_map, 'InitialMagnification', 'fit');

% Prikazivanje osa
axis on;


% Podesite oznake na osama (tick marks) i etikete (tick labels)
xticks(1:length(p_values));
xticklabels(arrayfun(@num2str, p_values, 'UniformOutput', false));
yticks(1:length(r_values));
yticklabels(arrayfun(@num2str, r_values, 'UniformOutput', false));

% Podešavanje oznaka i naslova grafika
xlabel('p');
ylabel('r');
title('HaLRTC');

saveas(gcf, 'graphRSE_HaLRTC.png');



