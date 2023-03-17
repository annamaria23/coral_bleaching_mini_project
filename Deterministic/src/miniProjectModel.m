function dY = miniProjectModel(t, Y, params)

% Parse parameters

n_row = params.dims(1);
n_col = params.dims(2);
i_max = params.dims(3);
j_max = params.dims(4);

p_HZ = params.H.p_HZ;
d_H = params.H.d_H;
a_H = params.H.a_H;

a_Z = params.Z.a_Z;
d_Z = params.Z.d_Z;
gamma = params.Z.gamma;
T_min = params.Z.T_min;

a_ROS = params.ROS.a_ROS;
D = params.ROS.D;
T_thres = params.ROS.T_thres;
L_thres = params.ROS.L_thres;
d_ROS_H = params.ROS.d_ROS_H;
d_ROS = params.ROS.d_ROS;

T_pattern = params.T.T_pattern;
T_initial = params.T.T_initial;
L_pattern = params.L;


% Initialise derivatives
n_grid = n_row * n_col;

dH = zeros(n_row, n_col);
dZ = zeros(n_row, n_col);
dROS = zeros(n_row, n_col);


% Get current function values
H = Y(1:n_grid);
Z = Y(n_grid + (1:n_grid));
ROS = Y(2*n_grid + (1:n_grid));

H = reshape(H, [n_row, n_col]);
Z = reshape(Z, [n_row, n_col]);
ROS = reshape(ROS, [n_row, n_col]);




% Calculate derivatives for the matrix (i, j) = (y, x)
for i = 1:n_row

    for j = 1:n_col


        T = get_temperature(t, T_pattern(i, j), T_initial);
        L = get_temperature(t, L_pattern(i, j), L_pattern(i, j));

        % ------- H equaitons ------

        dH(i,j) = a_H * p_HZ * Z(i,j) * H(i,j) * (1 - H(i,j)) - d_H * H(i,j);
        %         if H(i,j) > 0.01
        %             dH(i,j) = p_HZ * Z(i,j) * H(i,j) * (1 - H(i,j)) - d_H * H(i,j);
        %         else
        %             % Boundaries
        %             if i == 1
        %                 if j == 1
        %                     dH(i,j) = a_H * (H(i+1,j) + H(i,j+1)) / 2;
        %                 elseif j == n_col
        %                     dH(i,j) = a_H * (H(i+1,j) + H(i,j-1)) / 2;
        %                 else
        %                     dH(i,j) = a_H * (H(i+1,j) + H(i,j-1) + H(i,j+1)) / 3;
        %                 end
        %
        %             elseif i == n_row
        %                 if j == 1
        %                     dH(i,j) = a_H * (H(i-1,j) + H(i,j+1)) / 2;
        %                 elseif j == n_col
        %                     dH(i,j) = a_H * (H(i-1,j) + H(i,j-1)) / 2;
        %                 else
        %                     dH(i,j) = a_H * (H(i-1,j) + H(i,j-1) + H(i,j+1)) / 3;
        %                 end
        %
        %             elseif j == 1
        %                 dH(i,j) = a_H * (H(i-1,j) + H(i+1,j) + H(i,j+1)) / 3;
        %
        %             elseif j == n_col
        %                 dH(i,j) = a_H * (H(i-1,j) + H(i+1,j) + H(i,j-1)) / 3;
        %                 % Not at boundaries
        %             else
        %                 dH(i,j) = a_H * (H(i-1,j) + H(i+1,j) + H(i,j-1) + H(i,j+1)) / 4;
        %             end
        %         end




        % ------ Z equaitons ------
        if H(i,j) ~= 0 && Z(i, j) == 0
            dZ(i,j) = a_Z * H(i,j);
        else
            dZ(i,j) = a_Z * L * (T - T_min) * (1 - p_HZ) * Z(i,j) * (H(i,j) - Z(i,j)) - d_Z * Z(i,j) - gamma * Z(i, j) * ROS(i,j);
        end


        % ------ ROS equations ------


        if i == 1
            if j == 1
                dROS(i,j) = D * (ROS(i+1,j) + ROS(i,j+1) - 3 * ROS(i,j));
            elseif j == n_col
                dROS(i,j) = D * (ROS(i+1,j) + ROS(i,j-1) - 3 * ROS(i,j));
            else
                dROS(i,j) = D * (ROS(i+1,j) + ROS(i,j-1) + ROS(i,j+1) - 4 * ROS(i,j));
            end

        elseif i == n_row
            if j == 1
                dROS(i,j) = D * (ROS(i-1,j) + ROS(i,j+1) - 3 * ROS(i,j));
            elseif j == n_col
                dROS(i,j) = D * (ROS(i-1,j) + ROS(i,j-1) - 3 * ROS(i,j));
            else
                dROS(i,j) = D * (ROS(i-1,j) + ROS(i,j-1) + ROS(i,j+1) - 4 * ROS(i,j));
            end

        elseif j == 1
            dROS(i,j) = D * (ROS(i-1,j) + ROS(i+1,j) + ROS(i,j+1) - 4 * ROS(i,j));

        elseif j == n_col
            dROS(i,j) = D * (ROS(i-1,j) + ROS(i+1,j) + ROS(i,j-1) - 4 * ROS(i,j));

        else
            dROS(i,j) = D * (ROS(i-1,j) + ROS(i+1,j) + ROS(i,j-1) + ROS(i,j+1) - 5 * ROS(i,j));
        end

        dROS(i,j) = dROS(i,j) + a_ROS * max((T - T_thres), 0) * max((L - L_thres), 0) * Z(i,j) - d_ROS_H * H(i,j) * ROS(i,j) - d_ROS * T * ROS(i,j);

%         dROS(i,j) = a_ROS * max((T - T_thres), 0) * max((L - L_thres), 0) * Z(i,j) - d_ROS_H * H(i,j) * ROS(i,j) - d_ROS * T * ROS(i,j);

    end


end

dY = [dH(:); dZ(:); dROS(:)];


end