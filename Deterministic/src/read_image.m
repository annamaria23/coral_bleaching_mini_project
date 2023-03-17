function [new_H_matrix, new_T_matrix, img_matrix] = read_image(image, T_scale_image, padding, T_max, T_min, threshold)


X = imread(image);
X_scale = imread(T_scale_image);
X_scale = X_scale(:, 1, :);


X_scale_mat = reshape(X_scale, [size(X_scale, 1), 3]);
T_map = linspace(T_max, T_min, size(X_scale, 1));

n_row = 28;
n_col = 30;



col_increment = floor(size(X, 2) / n_col);
row_increment = floor(size(X, 1) / n_row);


img_matrix = zeros(n_row, n_col, 3);
T_matrix = zeros(n_row, n_col);
H_matrix = zeros(n_row, n_col);

new_H_matrix = zeros(n_row + 2*padding, n_col + 2*padding);
new_T_matrix = zeros(n_row + 2*padding, n_col + 2*padding);


for i = 1:n_row
    for j = 1:n_col

        img_matrix(i, j, :) = mean( ...
            X( ...
            (1:row_increment) + (i-1)*row_increment, ...
            (1:col_increment) + (j-1)*col_increment, ...
            : ...
            ), [1, 2]) ./ 255;

        %         if reshape(img_matrix(i, j, :), [1 3]) == [255 255 255]
        %             T_matrix(i, j) = NaN;
        %             continue
        %         end

        distance = reshape(img_matrix(i, j, :), [1, 3]) - double(X_scale_mat) ./ 255;
        distance = sum(distance.^2, 2);

        [~, idx] = min(distance);

        if img_matrix(i, j, 1) < threshold
            T_matrix(i, j) = 0;
            H_matrix(i, j) = 0;
        else
            T_matrix(i, j) = T_map(idx);
            H_matrix(i, j) = 1;

        end

        new_T_matrix(i + padding, j + padding) = T_matrix(i, j);
        new_H_matrix(i + padding, j + padding) = H_matrix(i, j);

    end
end

% imshow(img_matrix)
% 
% figure
% h = heatmap(new_T_matrix);
% h.ColorLimits = [-40 0];
% title('T pattern')
% 
% figure
% heatmap(new_H_matrix)
% title('H initial')






end
