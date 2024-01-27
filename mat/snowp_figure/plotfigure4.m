clear all;
clc
mat_dir = addPaths;
nrn_model_ver = 'maxH_Default';
cell_labels = string(["L1"; "L2/3"; "L4"; "L5"; "L6"]);
data_fold = fullfile(mat_dir, "snowp_data/1mV_10Hz");
phi = load(fullfile(data_fold, "phi.mat")); phi = phi.phi;
theta = load(fullfile(data_fold, "theta.mat")); theta = theta.theta;
subcell_peak_polar = cell(25, 1);

for cell_id=1:25
    data_file = sprintf('polarization_cell%g.mat', cell_id);
    if ~exist(fullfile(data_fold, data_file), "file")
        subpeakpolar = zeros(length(subcellular_index), 1);
        subcell_peak_polar{cell_id, 1} = subpeakpolar;
        continue
    end
    data = load(fullfile(data_fold, data_file));
    polarize266 = data.polarize266;

    secnames_fold = fullfile(mat_dir, "cell_data", nrn_model_ver, "secnames");
    secnames = load(fullfile(secnames_fold, sprintf('secnames%g.mat', cell_id)));
    secnames = secnames.secnames;
    secname = cellfun(@(x) split(x, '.'), secnames, 'UniformOutput', false);
    secname = cellfun(@(x) x{end-2}, secname, 'UniformOutput', false);
    secname = unique(secname, 'stable'); % senction name

    % soma, axon, Node, Myelin, Unmyelin, dend, apic
    subcellular_names = ["soma" "axon" "Node" "Myelin" "Unmyelin" "dend" "apic"]';
    subcellular_index = cell(length(subcellular_names), 1);
    for i=1:length(subcellular_index)
        temp_index = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(i)));
        subcellular_index{i} = temp_index;
    end

    polar_soma = cellfun(@(x) x(subcellular_index{1}), polarize266);
    [~, EF_direction_peaksomapolar] = max(polar_soma);

    subpeakpolar = zeros(length(subcellular_index), 1);
    for m=1:length(subcellular_index)
        if sum(subcellular_index{m}) ~= 0
            temppolar = cellfun(@(x) x(subcellular_index{m}), polarize266, 'UniformOutput', false);
            maxpolar1 = cellfun(@(x) max(x), temppolar, 'UniformOutput', false);
            subpeakpolar(m, 1) = max(cell2mat(maxpolar1));
        end
    end
    subcell_peak_polar{cell_id, 1} = subpeakpolar;
end


tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'soma'));
peakpolar_soma = cellfun(@(x) x(logical(tempindex)), subcell_peak_polar);
peakpolar_soma = reshape(peakpolar_soma, 5, 5);

tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'axon'));
tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'Node')) + tempindex; 
tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'Myelin')) + tempindex; 
tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'Unmyelin')) + tempindex; 
peakpolar_axon = cellfun(@(x) max(x(logical(tempindex))), subcell_peak_polar);
peakpolar_axon = reshape(peakpolar_axon, 5, 5);

tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'dend'));
peakpolar_basal = cellfun(@(x) x(logical(tempindex)), subcell_peak_polar);
peakpolar_basal = reshape(peakpolar_basal, 5, 5);

tempindex = cellfun(@(x) ~isempty(x), strfind(cellstr(subcellular_names), 'apic'));
peakpolar_apic = cellfun(@(x) x(logical(tempindex)), subcell_peak_polar);
peakpolar_apic = reshape(peakpolar_apic, 5, 5);

markers = {'+','^','x','o','s'}; 
colors = jet(5);

figdata = {peakpolar_soma; peakpolar_axon; peakpolar_basal; peakpolar_apic};
%% save data for python plotfigures
data_figure_peakpolar = struct();
data_figure_peakpolar.peakpolar_soma = peakpolar_soma;
data_figure_peakpolar.peakpolar_axon = peakpolar_axon;
data_figure_peakpolar.peakpolar_basal = peakpolar_basal;
data_figure_peakpolar.peakpolar_apic = peakpolar_apic;
% savefile = fullfile(mat_dir, 'snowp_data', 'Data_Figure4.mat');
% save(savefile, "data_figure_peakpolar");

%% plot figure4
t = tiledlayout(1,4);
for m=1:4
    nexttile
    data = figdata{m, 1};
    for i=1:5
        for j=1:5
            if data(j, i) ~= 0
                scatter(i, data(j, i), 60, Marker=markers(j),ColorVariable=colors(i), LineWidth=1.5);
                hold on;
            end
        end
    end
    ax = gca;
    ax.XLim = [0.5, 5.5];
    ax.YLim = [0, inf];
    ax.YGrid = 'on';
    ax.Box = 'off';
    ax.XTick = 1:5;
    ax.XTickLabel = cell_labels;
    ylabel('Peak polarization (mV)'); 
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

