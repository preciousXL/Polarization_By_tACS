function [secindex_soma, secindex_axon, secindex_basal, secindex_apic] = find_subcell_index(cell_id, mat_dir, nrn_model_ver)
    if nargin < 2
        mat_dir = addPaths;
        nrn_model_ver = 'maxH_Default';
    end
    secnames_fold = fullfile(mat_dir, 'cell_data', nrn_model_ver, 'secnames');
    secnames_file = fullfile(secnames_fold, sprintf('secnames%g.mat', cell_id));
    secnames = load(secnames_file);
    secnames = secnames.secnames;
    secname = cellfun(@(x) split(x, '.'), secnames, 'UniformOutput', false);
    secname = cellfun(@(x) x{end-2}, secname, 'UniformOutput', false);
    secname = unique(secname, 'stable');
    subcellular_names = ["soma" "axon" "Node" "Myelin" "Unmyelin" "dend" "apic"]';
    secindex_soma = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(1)));
    secindex_axon = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(2)));
    secindex_axon = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(3))) + secindex_axon;
    secindex_axon = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(4))) + secindex_axon;
    secindex_axon = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(5))) + secindex_axon;
    secindex_axon = logical(secindex_axon);
    secindex_basal = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(6)));
    secindex_apic = cellfun(@(x) ~isempty(x),strfind(secname,subcellular_names(7)));

end