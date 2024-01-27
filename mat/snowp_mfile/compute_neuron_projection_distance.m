% clearvars; clc;
%% 定义文件位置
mat_dir = 'H:\2Aberra_tACS\mat';
coord_fold = fullfile(mat_dir, 'cell_data/maxH_default/coordinates');
secname_fold = fullfile(mat_dir, 'cell_data/maxH_default/secnames');
%% 定义投影方向
list_phi = 0:15:360-1;
list_theta = 0:15:180;
list_direction = [list_theta(1) 0];
for i=2:length(list_theta)-1
    for j=1:length(list_phi)
        list_direction = [list_direction; [list_theta(i) list_phi(j)]];
    end
end
list_direction = [list_direction; [list_theta(end) 0]];

dist_proj_cell = cell(25, 1);
%% 求神经元每个segment在指定方向上的投影距离（soma为坐标原点[0, 0, 0]）
for i=1:25
    cell_id = i;
    % load coordinates and secnames of each neuron
    coord_filename = fullfile(coord_fold, sprintf('coordinates%g.mat', cell_id));
    secname_filename = fullfile(secname_fold, sprintf('secnames%g.mat', cell_id));
    coordnates = load(coord_filename);
    coordnates = coordnates.C;
    % data = [max_proj, min_proj, 有效计划长度L=|max_proj|+|min_proj|, 非对称性abs(|max_proj|-|min_proj|)/L]
    dist_proj = zeros(length(list_direction), 4);
    for j=1:length(list_direction)
        % Efield direction
        theta = list_direction(j, 1);
        phi = list_direction(j, 2);
        theta = deg2rad(theta);
        phi = deg2rad(phi);
        % compute each section projection on Efield direction
        unit_vec = [sin(theta)*cos(phi); 
                    sin(theta)*sin(phi); 
                    cos(theta)];
        seg_dot = coordnates*unit_vec;
        % seg_2norm = (sum(coordnates.^2, 2)).^0.5;
        seg_proj = seg_dot/1; % 在单位向量上的投影
        % data save
        max_proj = max(seg_proj);
        min_proj = min(seg_proj);
        Lpolar_proj = abs(max_proj) + abs(min_proj);
        asym_proj = abs(abs(max_proj) - abs(min_proj))/Lpolar_proj;
        dist_proj(j, :) = [max(abs(max_proj),abs(min_proj)), min(abs(max_proj),abs(min_proj)), Lpolar_proj, asym_proj];
    end
    dist_proj_cell{cell_id, 1} = dist_proj;
end
% data save
theta = list_direction(:, 1);
phi = list_direction(:, 2);
savefile = fullfile(mat_dir, 'snowp_data', 'dist_proj_25cells_default.mat');
save(savefile, "dist_proj_cell","theta","phi");


