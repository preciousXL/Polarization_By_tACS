function Edirect = find_maxpolar_Edirect(cell_id, subcellname, datafold)
    % cell_id: 1-25
    % subcellname: "soma" "axon" "Node" "Myelin" "Unmyelin" "dend" "apic"
    mat_dir = addPaths;
    if nargin < 3
        datafold = [mat_dir filesep 'snowp_data' filesep '1mV_10Hz'];
    end
    secnames_fold = fullfile(mat_dir, "cell_data", "maxH", "secnames");
    secnames = load(fullfile(secnames_fold, sprintf('secnames%g', cell_id)));
    secnames = secnames.secnames; 
    secname = cellfun(@(x) split(x, '.'), secnames, 'UniformOutput', false);
    secname = cellfun(@(x) x{end-2}, secname, 'UniformOutput', false);
    secname = unique(secname, 'stable'); % senction name
    subcell_index = cellfun(@(x) ~isempty(x),strfind(secname,subcellname));

    phi = load(fullfile(datafold, "phi.mat")); phi = phi.phi;
    theta = load(fullfile(datafold, "theta.mat")); theta = theta.theta;
    datafile = fullfile(datafold, sprintf('polarization_cell%g.mat', cell_id));
    data = load(datafile);
    polarize266 = data.polarize266;
    polarvalue = cellfun(@(x) x(subcell_index), polarize266, 'UniformOutput', false);
    maxpolar = cellfun(@(x) max(x), polarvalue);
    [~, Edirect_index] = max(maxpolar);

    Edirect = [theta(Edirect_index); phi(Edirect_index)];
end