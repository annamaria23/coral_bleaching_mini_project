


    % --------- Specify parameters -------------

    [~, T_pattern, ~] = read_image('Coral_temperature.png', 'Temperature_scale.png', 2, 31, 29, 0.7);
    [~, H_pattern, ~] = read_image('Coral_init.png', 'Coral_scale.png', 2, 0.8, 0, 0.25);
    [~, L_pattern, ~] =  read_image('Coral_depth.png', 'Depth_scale.png', 2, 0, -40, 0.2);

    L_pattern = (L_pattern + 40 * (L_pattern < 0)) ./ 40;


    effective_tiles = T_pattern > 0 & H_pattern > 0 & L_pattern > 0;
    n_effective = sum(effective_tiles, 'all');
    L_pattern(~effective_tiles) = 0;
    H_pattern(~effective_tiles) = 0;
    T_pattern(~effective_tiles) = 0;


    n_row = size(H_pattern, 1);
    n_col = size(H_pattern, 2);




    T_initial = 26;



    i_max = 1;
    j_max = 1;

    t_start = 0;
    t_end = 50;




    H_params.p_HZ = 0.4;  % Nutrient translocation
    H_params.d_H = 0.02;  % Host death rate
    H_params.a_H = 0.35;  % Host growth rate

    Z_params.a_Z = 0.11;  % Zoox growth rate
    Z_params.d_Z = 0.01;  % Zoox death rate
    Z_params.gamma = 250;  % ROS expell rate
    Z_params.T_min = 20;

    ROS_params.a_ROS = 1;  % ROS production rate

    ROS_params.D = 80;  % ROS diffusion constant (mm^2/day)
    ROS_params.T_thres = 27;  % Temperature threshold for ROS production
    ROS_params.L_thres = 0.5;  % Light threshold for ROS production
    ROS_params.d_ROS_H = 0.001;
    ROS_params.d_ROS = 0.005;





    % Loop through each proportion (can be only one proportion)

    disp('==========================================')

    % --------- Params to the function ---------

    params.dims = [n_row, n_col, i_max, j_max];
    params.H = H_params;
    params.Z = Z_params;
    params.ROS = ROS_params;
    params.T.T_pattern = T_pattern;
    params.T.T_initial = T_initial;
    params.L =  L_pattern;



    % --------- Initialisation ---------
    n_grid = n_row * n_col;
    H_initial = H_pattern(:);
    initial = [H_initial; H_initial; zeros(n_grid, 1)];


    % --------- Solve the ODE ---------
    options = odeset('RelTol', 1e-4, 'AbsTol',1e-6);
    [t, Y_solution] = ode45(@(t, Y) miniProjectModel(t, Y, params), t_start:t_end, initial, options);

    disp('Solved!!')
    disp('==========================================')



    % --------- Interpret ODE solution ---------

    H_vec = Y_solution(:, 1:n_grid);
    Z_vec = Y_solution(:, n_grid + (1:n_grid));
    ROS_vec = Y_solution(:, 2*n_grid + (1:n_grid));

    H = zeros(n_row, n_col, length(t));
    Z = zeros(n_row, n_col, length(t));
    ROS = zeros(n_row, n_col, length(t));
    T_mat = zeros(n_row, n_col, length(t));
    L_mat = zeros(n_row, n_col, length(t));


    for i = 1:n_row
        for j = 1:n_col
            for k = 1:length(t)
                T_mat(i, j, k) = get_temperature(t(k), params.T.T_pattern(i, j), params.T.T_initial) * effective_tiles(i,j);
                L_mat(i, j, k) = get_temperature(t(k), params.L(i, j), params.L(i, j)) * effective_tiles(i,j);
            end
        end
    end



    H_sum = zeros(length(t), 1);
    Z_sum = zeros(length(t), 1);
    ROS_sum = zeros(length(t), 1);
    T = zeros(length(t), 1);

    for i = 1:length(t)

        H(:, :, i) = reshape(H_vec(i, :), [n_row, n_col]);
        Z(:, :, i) = reshape(Z_vec(i, :), [n_row, n_col]);
        ROS(:, :, i) = reshape(ROS_vec(i, :), [n_row, n_col]);

        H_sum(i, 1) = sum(H(:, :, i), 'all');
        Z_sum(i, 1) = sum(Z(:, :, i), 'all');
        ROS_sum(i, 1) = sum(ROS(:, :, i), 'all');

        T_mean_mat = T_mat(:, :, i);

        T(i, 1) = mean(T_mean_mat(effective_tiles), 'all');


    end


    Solution.H = H;
    Solution.Z = Z;
    Solution.ROS = ROS;
    Solution.T_mat = T_mat;
    Solution.L_mat = L_mat;



    % ----------------- Plot --------------

    fig = figure;

    subplot(2,1,2)
    yyaxis left
    plot(t, H_sum/n_effective, t, Z_sum/n_effective, 'LineWidth', 2)
    ylabel('Host and Zoox')

    yyaxis right
    plot(t, ROS_sum/n_effective, 'LineWidth', 2)
    ylabel('ROS')


    xlabel('t')
    title('System')
    legend('Host', 'Zoox', '[ROS]')
    set(gca, 'FontWeight', 'bold', 'FontSize', 12)


    subplot(2,1,1)
    plot(t, T, 'LineWidth', 2)
    yline(ROS_params.T_thres, 'r-', 'LineWidth', 1, 'Label', 'Threshold')
    xlabel('t')
    ylabel('Temperature')
    title('Temperature')
    ylim([0, 40])
    set(gca, 'FontWeight', 'bold', 'FontSize', 12)

    sgtitle('L-thres = ' + string(params.ROS.L_thres))







    
%     save_path = fullfile(pwd, 'Results', string(datetime) + '- p' + string(p));
%     save_path = fullfile(pwd, 'Results', string(datetime));
% 
%     
%     
%     if ~exist(save_path, "dir")
%         mkdir(save_path);
%     end
%     
%     addpath(fullfile(pwd, 'Results'))
%     addpath(save_path);
%     
%     
%     
%     savefig(fig, fullfile(save_path, 'LinePlot'), "compact")
%     save(fullfile(save_path, 'Params'), 'params')
%     save(fullfile(save_path, 'Solution'), "Solution")
%     save(fullfile(save_path, 'Initial'), "initial")
%     
%     
%     plot_matrix = [t, H_sum/n_effective, Z_sum/n_effective, ROS_sum/n_effective];
%     writematrix(plot_matrix, fullfile(save_path, 'LinePlot_Data.csv'))




%%
% Arraydata = table2array(Arraydata);
% p_array = Arraydata(:, 1);
% H_array = Arraydata(:, 2);
% Z_array = Arraydata(:, 3);
% H_p_array = Arraydata(:, 4);
% Z_p_array = Arraydata(:, 5);

figure
plot(p_array, H_array/400, p_array, Z_array/400, p_array, H_p_array, p_array, Z_p_array, 'LineWidth', 2)
legend('Coral', 'Algae', 'Coral per unit', 'Algae per unit')

% Array_data = [p_array', H_array', Z_array', H_p_array', Z_p_array'];
% writematrix(Array_data, fullfile(save_path, 'Array_data.csv'))

xlabel('Initial coral density')
ylabel('End time point density')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')

%% ---------------- Animation --------------------
set(0, 'DefaultFigureRenderer', 'painters')


% figure
%     set(gca,'FontSize',12,'FontWeight','bold')



H = Solution.H .* 100;
Z = Solution.Z .* 100;
ROS = Solution.ROS;
T_mat = Solution.T_mat;
L_mat = Solution.L_mat;
t = 1:size(H, 3);





c = jet(100);
c_map = c(10:90, :);

% for i = 30 : length(t)
for i = 1:length(t)



    % H
    current_H = H(:, : ,i);
    current_H(~effective_tiles) = NaN;

    subplot(1, 5, 4)
    h = heatmap(current_H, 'CellLabelColor', 'auto', 'CellLabelFormat', '%.2f');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));

    if min(H, [], "all") == max(H, [], "all")
        h.ColorLimits = [0, 1];
    elseif true
        h.ColorLimits = [0, 100];
        h.Colormap = c_map;
    elseif min(H, [], "all") < 0 || max(H, [], "all") > 0
        h.ColorLimits = [min(H, [], "all") max(H, [], "all")];

    end
    title('Coral (%)')




    % Z
    current_Z = Z(:, : ,i);
    current_Z(~effective_tiles) = NaN;

    subplot(1, 5, 5)
    h = heatmap(current_Z, 'CellLabelColor', 'none');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    h.ColorLimits = [0, 1];
    if min(Z, [], "all") == max(Z, [], "all")
        h.ColorLimits = [0, 1];
    elseif true
        h.ColorLimits = [0, 100];
        h.Colormap = c_map;

    elseif min(Z, [], "all") < 0 || max(Z, [], "all") > 1
        h.ColorLimits = [min(Z, [], "all") max(Z, [], "all")];
    else
        h.ColorLimits = [0, 1];
    end
    title('Zooxanthellae (%)')




    % ROS
    current_ROS = ROS(:, : ,i);
    current_ROS(~effective_tiles) = NaN;

    subplot(1, 5, 3)
    h = heatmap(current_ROS, 'CellLabelColor', 'none');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    if min(ROS, [], "all") == max(ROS, [], "all")
        h.ColorLimits = [0, 1];
    elseif true
        h.ColorLimits = [0.001 0.005];
        h.ColorLimits = [min(ROS, [], "all") max(ROS, [], "all")];
        h.Colormap = c_map;

    else
    end
    title('ROS (a.u.)')



    % T_mat
    current_T_mat = T_mat(:, : ,i);
    current_T_mat(~effective_tiles) = NaN;

    subplot(1, 5, 1)
    h = heatmap(current_T_mat, 'CellLabelColor', 'none');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    if min(T_mat, [], "all") == max(T_mat, [], "all")
        h.ColorLimits = [0, 1];
    elseif true
        h.ColorLimits = [27 31];
        h.Colormap = c_map;
    else
        h.ColorLimits = [min(T_mat, [], "all") max(T_mat, [], "all")];
    end
    title('Temperature (˚C)')


    % L_mat
    current_L_mat = L_mat(:, : ,i);
    current_L_mat(~effective_tiles) = NaN;

    subplot(1, 5, 2)
    h = heatmap(current_L_mat, 'CellLabelColor', 'none');
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    if min(L_mat, [], "all") == max(L_mat, [], "all")
        h.ColorLimits = [0, 1];
    elseif true
        h.ColorLimits = [0 1];
        h.Colormap = c_map;

    else
        h.ColorLimits = [min(T_mat, [], "all") max(T_mat, [], "all")];
    end
    title('Light Intensity (a.u.)')




    sgtitle('t = ' + string(i))

    drawnow

end







