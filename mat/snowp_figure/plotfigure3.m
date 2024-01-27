mat_dir = addPaths;
order = 1:1:25;
order = reshape(order, 5, 5);
order = reshape(order', 25, 1);

data_figure_maxVmDeflection = cell(25, 1);

figure('Units','centimeters','Position',[10, 10, 15, 10])
t = tiledlayout(5,5);
for i=1:25
    figdata = struct();
    nexttile
    cell_id = order(i);
    datafold = fullfile(mat_dir, 'snowp_data', '10Hz_amp');
    datafile = fullfile(datafold, sprintf('polarization_cell%g.mat', cell_id));
    if ~exist(datafile, "file")
        continue
    end
    data = load(datafile);
    list_amp = data.list_amp;
    polarize_amp = data.maxVmshift; 

    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_index(cell_id, mat_dir, 'maxH_default');
    colors = ['k', 'r', 'b', 'g'];
    polar_soma  = cellfun(@(x) x(idxsoma), polarize_amp, 'UniformOutput', false);
    polar_axon  = cellfun(@(x) x(idxaxon), polarize_amp, 'UniformOutput', false);
    polar_basal = cellfun(@(x) x(idxbasal), polarize_amp, 'UniformOutput', false);
    polar_apic  = cellfun(@(x) x(idxapic), polarize_amp, 'UniformOutput', false); 
    [~, idxsoma_sec] = cellfun(@(x) max(x), polar_soma); 
    [~, idxaxon_sec] = cellfun(@(x) max(x), polar_axon); 
    [~, idxbasal_sec] = cellfun(@(x) max(x), polar_basal); 
    idxsoma_sec = mode(idxsoma_sec);
    idxaxon_sec = mode(idxaxon_sec);
    idxbasal_sec = mode(idxbasal_sec);

    polar_soma = cellfun(@(x) x(idxsoma_sec), polar_soma);
    polar_axon = cellfun(@(x) x(idxaxon_sec), polar_axon);
    polar_basal = cellfun(@(x) x(idxbasal_sec), polar_basal);

    figdata.polar_soma = polar_soma;
    figdata.polar_axon = polar_axon;
    figdata.polar_basal = polar_basal;
    flagnorm = 0;
    if ~isempty(polar_apic{1})
        [~, idxapic_sec] = cellfun(@(x) max(x), polar_apic); 
        idxapic_sec = mode(idxapic_sec);
        polar_apic = cellfun(@(x) x(idxapic_sec), polar_apic);
        figdata.polar_apic = polar_apic;
        if flagnorm == 1
            [polar_apic, ~] = mapminmax(polar_apic', 0, 1);
            polar_apic = polar_apic';
        end
        plot(list_amp, polar_apic, colors(4), LineWidth=2); hold on;
    end
    if flagnorm == 1
        [polar_soma, ~] = mapminmax(polar_soma', 0, 1);
        polar_soma = polar_soma';
        [polar_axon, ~] = mapminmax(polar_axon', 0, 1);
        polar_axon = polar_axon';
        [polar_basal, ~] = mapminmax(polar_basal', 0, 1);
        polar_basal = polar_basal';
    end
    plot(list_amp, polar_soma, colors(1), LineWidth=2); hold on;
    plot(list_amp, polar_axon, colors(2), LineWidth=2); hold on;
    plot(list_amp, polar_basal, colors(3), LineWidth=2); hold on;
    xlim([0, 3]); 
    ylim([0, 3]);
    ax = gca;
    ax.XTick = [0, 1, 2, 3];
    ax.Box = 'off';
    if ismember(i, 1:20)
        ax.XTick = [];
    else
        ax.XLabel.String = 'Amplitude (mV/mm)';
    end
    if ismember(i, 1:5:25)
        ax.YLabel.String = 'Lp (mm)';
    end
    data_figure_maxVmDeflection{i, 1} = figdata;
end
t.TileSpacing = 'compact';
t.Padding = 'compact';


%% save data for python plotfigures
savefile = fullfile(mat_dir, 'snowp_data', 'Data_Figure_intensity.mat');
save(savefile, "data_figure_maxVmDeflection", "list_amp");