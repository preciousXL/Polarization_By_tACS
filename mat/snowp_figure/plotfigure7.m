%% all 25 cells
polarsoma25cells = [];
distproj25cells = [];
asymm25cells = [];

mat_dir = addPaths;
order = 1:1:25;
order = reshape(order, 5, 5);
order = reshape(order', 25, 1);

data_figure_sub_peakpolar = cell(25, 1);

t = tiledlayout(5,5);
for m=1:25
    cell_id = order(m);
    datafold_polarization = fullfile(mat_dir, "snowp_data/1mV_10Hz");
    datafile_polarization = fullfile(datafold_polarization, sprintf('polarization_cell%g.mat', cell_id));
    datafile_distproj = fullfile(mat_dir, "snowp_data/dist_proj_25cells_Default.mat");
    data = load(datafile_polarization);

    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_index(cell_id, mat_dir, 'maxH_Default');
    
    polarize266 = data.polarize266;
    polar_soma = zeros(length(polarize266), 1);
    polar_axon = zeros(length(polarize266), 1);
    polar_basal = zeros(length(polarize266), 1);
    polar_apic = zeros(length(polarize266), 1);
    polar_allsub = zeros(length(polarize266), 1);
    for i=1:length(polarize266)
        polar_soma(i, 1) = max(polarize266{i}(idxsoma)); 
        polar_axon(i, 1) = max(polarize266{i}(idxaxon)); 
        polar_basal(i, 1) = max(polarize266{i}(idxbasal)); 
        if sum(idxapic)~=0
            polar_apic(i, 1) = max(polarize266{i}(idxapic)); 
        else
            polar_apic(i, 1) = 0;
        end
        polar_allsub(i, 1) = max(polarize266{i}); 
    end
    subcellular = struct();
    subcellular.polar_soma = polar_soma;
    subcellular.polar_axon = polar_axon;
    subcellular.polar_basal = polar_basal;
    subcellular.polar_apic = polar_apic;
    subcellular.polar_allsub = polar_allsub;
    
    data = load(datafile_distproj);
    dist_proj_25cells = data.dist_proj_cell;
    dist_proj = dist_proj_25cells{cell_id, 1}(:, 3);
    cell_asymm = dist_proj_25cells{cell_id, 1}(:, 4);

    subcellular.dist_proj = dist_proj;
    subcellular.cell_asymm = cell_asymm;
    subcellular.maxproj = dist_proj_25cells{cell_id, 1}(:, 1);
    subcellular.minproj = dist_proj_25cells{cell_id, 1}(:, 2);
    data_figure_sub_peakpolar{m, 1} = subcellular;

    polarsoma25cells = [polarsoma25cells; polar_soma];
    distproj25cells = [distproj25cells; dist_proj];
    asymm25cells = [asymm25cells; cell_asymm];

    coeffs = polyfit(dist_proj, polar_soma, 1);
    tempx = linspace(min(dist_proj), max(dist_proj), 500);
    yfit = polyval(coeffs, dist_proj);
    R2 = 1 - sum((yfit - polar_soma).^2)/((length(yfit)-1)*var(polar_soma));
    nexttile
    scatter(dist_proj, polar_soma, 25, 'filled'); 
    hold on;
    plot(tempx, polyval(coeffs, tempx), 'r', LineWidth=2);
    title(['R2=' num2str(R2)])
end
t.TileSpacing = 'compact';
t.Padding = 'compact';

%% save data for python plotfigures
savefile = fullfile(mat_dir, 'snowp_data', 'Data_Figure_moreffect.mat');
save(savefile, "data_figure_sub_peakpolar");

