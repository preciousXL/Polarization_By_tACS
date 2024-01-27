%% phase-frequency characteristic
mat_dir = addPaths;
filepath = [mat_dir, filesep, 'snowp_data\maxpolar_index_1mV_10Hz.mat'];
load(filepath)

biophysics = '_Default';  % '_Exc', '_Inh'
orders = 1:1:25;
orders = reshape(orders, 5, 5);
orders = reshape(orders', 25, 1);
data_figure_phasefreq = cell(25, 1);
t = tiledlayout(5,5);
for i=1:25
    nexttile(i)
    cell_id = orders(i);
    datafold = fullfile(mat_dir, 'snowp_data', '1mV_freq', ['phase' biophysics]);
    datafile = fullfile(datafold, sprintf('polarization_cell%g.mat', cell_id));
    if ~exist(datafile, "file")
        continue
    end
    data = load(datafile);
    list_f = data.list_f;
    idxsoma_sec = maxpolar_index_1mV_10Hz(cell_id, 1);
    idxaxon_sec = maxpolar_index_1mV_10Hz(cell_id, 2);
    idxbasal_sec = maxpolar_index_1mV_10Hz(cell_id, 3);
    idxapic_sec = maxpolar_index_1mV_10Hz(cell_id, 4);

    polarize_phase = data.polarize_phase;
    colors = ['k', 'r', 'b', 'g'];
    polar_soma = cellfun(@(x) x(idxsoma_sec), polarize_phase);
    polar_axon = cellfun(@(x) x(idxaxon_sec), polarize_phase);
    polar_basal = cellfun(@(x) x(idxbasal_sec), polarize_phase);

    polar_soma = mod(polar_soma, 180);
    polar_axon = mod(polar_axon, 180);
    polar_basal = mod(polar_basal, 180);

    datapython = struct();
    datapython.phase_soma = polar_soma;
    datapython.phase_axon = polar_axon;
    datapython.phase_basal = polar_basal;
    % plot
    if idxapic_sec ~= 0
        polar_apic = cellfun(@(x) x(idxapic_sec), polarize_phase);
        polar_apic = mod(polar_apic, 180);
        datapython.phase_apic = polar_apic;
        plot(list_f, polar_apic, colors(4), LineWidth=2); hold on;
    end
    plot(list_f, polar_soma, colors(1), LineWidth=2); hold on;
    plot(list_f, polar_axon, colors(2), LineWidth=2); hold on;
    plot(list_f, polar_basal, colors(3), LineWidth=2); hold on;
    xlim([0, 100]); 
    % ylim([-220, 220]);
    data_figure_phasefreq{i, 1} = datapython;
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% save data for python plotfigures
savefile = fullfile(mat_dir, 'snowp_data', sprintf('Data_Figure_PhaseFreq%s.mat', biophysics));
save(savefile, "data_figure_phasefreq", "list_f");
