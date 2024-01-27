%% cell morphology，2D or 3D

data_mor = cell(25, 1);

mat_dir = addPaths;
nrn_dir = [mat_dir '/../nrn'];
nrn_model_ver = 'maxH'; 
iontype = 'Default';
nrn_model_ver = sprintf('%s_%s', nrn_model_ver, iontype);

order = 1:1:25;
order = reshape(order, 5, 5);
order = reshape(order', 25, 1);
figure('Units','centimeters','Position',[10, 10, 15, 10])
t = tiledlayout(5,5);
for i=1:25
    data = struct(); 
    nexttile(i);
    cell_id = order(i);
    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_comp(cell_id, mat_dir, nrn_model_ver);
    colors = ['k', 'r', 'b', 'g'];
   
    coordfold = fullfile(mat_dir, "cell_data", nrn_model_ver, "coordinates");
    parentfold = fullfile(mat_dir, "cell_data", nrn_model_ver, "parent_inds");
    coordfile = fullfile(coordfold, sprintf('coordinates%g.mat', cell_id));
    parentfile = fullfile(parentfold, sprintf('parent_inds%g.mat', cell_id));
    coordinates = load(coordfile); coordinates = coordinates.C;
    parent_inds = load(parentfile); parent_inds = parent_inds.parent_inds;
    
    data.coord_soma = coordinates(1, :);
    % axon
    num_comp = sum(idxaxon);
    morph_axon = zeros(3, num_comp*3);
    morph_axon(:, 1:3:end) = coordinates(parent_inds(idxaxon(2:end)), :)';
    morph_axon(:, 2:3:end) = coordinates(idxaxon, :)';
    morph_axon(:, 3:3:end) = nan(3, num_comp);
    data.coord_axon = morph_axon;
    % basal dend
    num_comp = sum(idxbasal);
    morph_basal = zeros(3, num_comp*3);
    morph_basal(:, 1:3:end) = coordinates(parent_inds(idxbasal(2:end)), :)';
    morph_basal(:, 2:3:end) = coordinates(idxbasal, :)';
    morph_basal(:, 3:3:end) = nan(3, num_comp);
    data.coord_basal = morph_basal;
    % apic dend
    if sum(idxapic) > 0
        num_comp = sum(idxapic);
        morph_apic = zeros(3, num_comp*3);
        morph_apic(:, 1:3:end) = coordinates(parent_inds(idxapic(2:end)), :)';
        morph_apic(:, 2:3:end) = coordinates(idxapic, :)';
        morph_apic(:, 3:3:end) = nan(3, num_comp);
        data.coord_apic = morph_apic;
    end
    % axon + apic + basal
    num_comp = size(coordinates, 1);
    morph_plot = zeros(3, (num_comp-1)*3); % shape=((x,y,z), numcomp*3-3)
    morph_plot(:, 1:3:end) = coordinates(parent_inds, :)';
    morph_plot(:, 2:3:end) = coordinates(2:end, :)';
    morph_plot(:, 3:3:end) = nan(3, num_comp-1);

    data_mor{i, 1} = data;
    % plot
    if false
        x = morph_plot(1, :);
        y = morph_plot(2, :);
        z = morph_plot(3, :);
        plot3(x, y, z, 'Color','grey'); hold on;
        scatter3(coordinates(1, 1), coordinates(1, 2), coordinates(1, 3), 25, 'k', 'filled')
        view(0, 0)
    else
        x = morph_axon(1, :); y = morph_axon(2, :); z = morph_axon(3, :); 
        hlineaxon = plot3(x, y, z, 'Color','red'); hold on;
        x = morph_basal(1, :); y = morph_basal(2, :); z = morph_basal(3, :); 
        hlinebasal = plot3(x, y, z, 'Color','green'); hold on;
        if sum(idxapic) > 0
            x = morph_apic(1, :); y = morph_apic(2, :); z = morph_apic(3, :); 
            hlineapic = plot3(x, y, z, 'Color','blue'); hold on;
        end
        hlinesoma = scatter3(coordinates(1, 1), coordinates(1, 2), coordinates(1, 3), 25, 'k', 'filled');
        view(0, 0)
    end
    if sum(idxapic) > 0
        set(gca, 'Children', [hlineaxon, hlinebasal, hlineapic, hlinesoma]);
    else
        set(gca, 'Children', [hlineaxon, hlinebasal, hlinesoma]);
    end
    ax = gca;
    ax.LineWidth = 1;
    
    % ax.Visible = 'off';
    axis tight
    axis equal
    if ismember(i, 1:5:25)
        ax.XLim = [-450, 300];
        ax.ZLim = [-220, 120];
    elseif ismember(i, 2:5:25)
        ax.XLim = [-1500, 1000];
        ax.ZLim = [-1300, 500];
    elseif ismember(i, 3:5:25)
        ax.XLim = [-900, 500];
        ax.ZLim = [-800, 700];
    elseif ismember(i, 4:5:25)
        ax.XLim = [-1000, 1200];
        ax.ZLim = [-1200, 1200];
    elseif ismember(i, 5:5:25)
        ax.XLim = [-800, 800];
        ax.ZLim = [-800, 1300];
    end
    ax.XTick = []; ax.ZTick = [];

end
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% save data for python plotfigures
savefile = fullfile(mat_dir, 'snowp_data', 'Data_Figure_morphology.mat');
save(savefile, "data_mor");



