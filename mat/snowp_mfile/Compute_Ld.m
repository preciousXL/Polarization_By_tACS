mat_dir = addPaths;
coord_fold = fullfile(mat_dir, 'cell_data/maxH_Default/coordinates');
secname_fold = fullfile(mat_dir, 'cell_data/maxH_Default/secnames');

list_phi = 0:15:360-1;
list_theta = 0:15:180;
list_direction = [list_theta(1) 0];
for i=2:length(list_theta)-1
    for j=1:length(list_phi)
        list_direction = [list_direction; [list_theta(i) list_phi(j)]];
    end
end
list_direction = [list_direction; [list_theta(end) 0]];

Ld_soma_axon_basal_apic_all = cell(25, 1); % 
for i=1:25
    cell_id = i;
    % load coordinates and secnames of each neuron
    coord_filename = fullfile(coord_fold, sprintf('coordinates%g.mat', cell_id));
    secname_filename = fullfile(secname_fold, sprintf('secnames%g.mat', cell_id));
    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_segment_index(cell_id);
    coordnates = load(coord_filename);
    coordnates = coordnates.C;
    dist_proj = zeros(length(list_direction), 5); 
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
        seg_proj = seg_dot/1; 
        % data save
        max_proj = max(seg_proj(idxsoma)); min_proj = min(seg_proj(idxsoma)); % soma
        dist_proj(j, 1) = abs(max_proj) + abs(min_proj);
        max_proj = max(seg_proj(idxaxon)); min_proj = min(seg_proj(idxaxon)); % axon
        dist_proj(j, 2) = abs(max_proj) + abs(min_proj);
        max_proj = max(seg_proj(idxbasal)); min_proj = min(seg_proj(idxbasal)); % basal
        dist_proj(j, 3) = abs(max_proj) + abs(min_proj);
        if sum(idxapic) ~= 0
            max_proj = max(seg_proj(idxapic)); min_proj = min(seg_proj(idxapic)); % apic
            dist_proj(j, 4) = abs(max_proj) + abs(min_proj);
        end
        max_proj = max(seg_proj); min_proj = min(seg_proj); % all
        dist_proj(j, 5) = abs(max_proj) + abs(min_proj);
    end
    Ld_soma_axon_basal_apic_all{cell_id, 1} = dist_proj;
end
% data save
theta = list_direction(:, 1);
phi = list_direction(:, 2);
savefile = fullfile(mat_dir, 'snowp_data', 'Revise_Data', 'Ld_soma_axon_basal_apic_all.mat');
save(savefile, "Ld_soma_axon_basal_apic_all", "theta", "phi");