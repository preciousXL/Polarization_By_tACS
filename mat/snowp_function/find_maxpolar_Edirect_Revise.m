function Edirect = find_maxpolar_Edirect_Revise(cell_id, subcellname, datafold)
    % cell_id: 1-25
    % subcellname: "soma" "axon" "apic"
    % axon: "axon" "Node" "Myelin" "Unmyelin" "dend"
    if nargin < 3
        datafold = 'H:\2Aberra_tACS\mat\snowp_data\1mV_10Hz';
    end
    mat_dir = 'H:\2Aberra_tACS\mat';
    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_index(cell_id, mat_dir, 'maxH_default');
    if subcellname == "soma"
        subcell_index = idxsoma;
    elseif subcellname == "axon"
        subcell_index = idxaxon;
    elseif subcellname == "basal"
        subcell_index = idxbasal;
    elseif subcellname == "apic"
        subcell_index = idxapic;
    end
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