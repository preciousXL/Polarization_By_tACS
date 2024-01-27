mat_dir = addPaths;
data_fold = fullfile(mat_dir, "snowp_data/1mV_10Hz");
phi = load(fullfile(data_fold, 'phi.mat')); phi = phi.phi;
theta = load(fullfile(data_fold, 'theta.mat')); theta = theta.theta;
list_theta = 0:15:180;
list_phi = 0:15:360-1;
iontype = 'Default';

order = 1:1:25;
order = reshape(order, 5, 5);
order = reshape(order', 25, 1);
t = tiledlayout(5,5);
for m=1:25
    cell_id = order(m);
    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_index(cell_id, mat_dir, sprintf('maxH_%s', iontype));
    if ~ismember(cell_id,  1:25)
        nexttile;
        continue;
    end
    data_file = sprintf('polarization_cell%g.mat', cell_id);
    data = load(fullfile(data_fold, data_file));
    polarize266 = data.polarize266;
    polarize266_phase = data.polarize_phase;
    
    polar_soma = zeros(length(polarize266), 1);
    polar_axon = zeros(length(polarize266), 1);
    polar_basal = zeros(length(polarize266), 1);
    polar_apic = zeros(length(polarize266), 1);
    for i=1:length(polarize266)
        polar_soma(i, 1) = polarize266{i}(idxsoma);
        if polarize266_phase{i}(idxsoma) >= 180
            polar_soma(i, 1) = - polar_soma(i, 1);
        end
        [polar_axon(i, 1), idx00] = max(polarize266{i}(idxaxon));
        if polarize266_phase{i}(idx00) >= 180
            polar_axon(i, 1) = - polar_axon(i, 1);
        end
        [polar_basal(i, 1), idx00] = max(polarize266{i}(idxbasal));
        if polarize266_phase{i}(idx00) >= 180
            polar_basal(i, 1) = - polar_basal(i, 1);
        end
        if sum(idxapic)~=0
            [polar_apic(i, 1), idx00] = max(polarize266{i}(idxapic));
            if polarize266_phase{i}(idx00) >= 180
                polar_apic(i, 1) = - polar_apic(i, 1);
            end
        end
    end

    % optimal field direction
    polar_soma = abs(polar_soma);
    polar_axon = abs(polar_axon);
    polar_basal = abs(polar_basal);
    [~, maxidx_soma] = max(polar_soma);
    [~, maxidx_axon] = max(polar_axon);
    [~, maxidx_basal] = max(polar_basal);
    disp([cell_id, theta(maxidx_basal)]);

    if sum(idxapic)~=0
        polar_apic = abs(polar_apic);
        [~, maxidx_apic] = max(polar_apic);
    end

    datadirect_soma = zeros(length(list_theta), length(list_phi));
    datadirect_axon = zeros(length(list_theta), length(list_phi));
    datadirect_basal = zeros(length(list_theta), length(list_phi));
    datadirect_apic = zeros(length(list_theta), length(list_phi));
    for ii=1:length(list_theta)
        for jj=1:length(list_phi)
            temptheta = list_theta(ii);
            if temptheta==0 || temptheta==180
                tempphi = 0;
            else
                tempphi = list_phi(jj);
            end
            idx_phi = phi==tempphi;
            idx_theta = theta==temptheta;
            tempidx = idx_phi & idx_theta;
            datadirect_soma(ii, jj) = polar_soma(tempidx);
            datadirect_axon(ii, jj) = polar_axon(tempidx);
            datadirect_basal(ii, jj) = polar_basal(tempidx);
            datadirect_apic(ii, jj) = polar_apic(tempidx); 
        end
    end
    % index for -90°
    idx_270 = find(0:15:360-1 == 90);
    % soma
    tempdata1 = datadirect_soma(:, idx_270:end);
    tempdata2 = datadirect_soma(:, 1:idx_270-1);
    tempdata = [tempdata1 tempdata2];
    tempdata = [tempdata tempdata(:, 1)];
    tempdata = interp2(tempdata, 4);
    if m==1
        data_figure_direction_soma = zeros(25, size(tempdata, 1), size(tempdata, 2));
        data_figure_direction_axon = zeros(25, size(tempdata, 1), size(tempdata, 2));
        data_figure_direction_basal = zeros(25, size(tempdata, 1), size(tempdata, 2));
        data_figure_direction_apic = zeros(25, size(tempdata, 1), size(tempdata, 2));
    end
    data_figure_direction_soma(m, :, :) = tempdata;
    tempdata = interp2(tempdata, 4);
    clear tempdata tempdata1 tempdata2
    % axon
    tempdata1 = datadirect_axon(:, idx_270:end);
    tempdata2 = datadirect_axon(:, 1:idx_270-1);
    tempdata = [tempdata1 tempdata2];
    tempdata = [tempdata tempdata(:, 1)];
    tempdata = interp2(tempdata, 4);
    data_figure_direction_axon(m, :, :) = tempdata;
    clear tempdata tempdata1 tempdata2
    % basal
    tempdata1 = datadirect_basal(:, idx_270:end);
    tempdata2 = datadirect_basal(:, 1:idx_270-1);
    tempdata = [tempdata1 tempdata2];
    tempdata = [tempdata tempdata(:, 1)];
    tempdata = interp2(tempdata, 4);
    data_figure_direction_basal(m, :, :) = tempdata;
    clear tempdata tempdata1 tempdata2
    % apic
    tempdata1 = datadirect_apic(:, idx_270:end);
    tempdata2 = datadirect_apic(:, 1:idx_270-1);
    tempdata = [tempdata1 tempdata2];
    tempdata = [tempdata tempdata(:, 1)];
    tempdata = interp2(tempdata, 4);
    data_figure_direction_apic(m, :, :) = tempdata;
    clear tempdata tempdata1 tempdata2
    
    % normalization
    % [soma_norm, ~] = mapminmax(polar_soma', 0, 1);
    [soma_norm, ~] = mapminmax(polar_soma', -1, 1);
    soma_norm = soma_norm';
    
    cmap = fake_parula(1000); 
    display_scale = 'lin';
    nexttile
    % fig = figure('units','normalized','outerposition',[0.2 0.2 0.4 0.35],'Color','w');
    plotThreshMapMW(soma_norm,theta,phi, 'center_phi', -90);
    ax = gca;
    colormap(ax,cmap); 
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% save data for python plotfigures
list_theta = deg2rad(linspace(-90, 90, size(squeeze(data_figure_direction_soma(1, :, :)), 1)));
list_phi = deg2rad(linspace(-180, 180, size(squeeze(data_figure_direction_soma(1, :, :)), 2)));
savefile = fullfile(mat_dir, 'snowp_data', 'Data_Figure_direction.mat');
save(savefile, "data_figure_direction_apic", "list_theta", "list_phi", "data_figure_direction_basal", "data_figure_direction_axon", "data_figure_direction_soma");



