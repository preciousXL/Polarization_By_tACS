function [segindex_soma, segindex_axon, segindex_basal, segindex_apic] = find_subcell_segment_index(cell_id, mat_dir, nrn_model_ver)
    if nargin < 2
        mat_dir = 'H:\2Aberra_tACS\mat';
        nrn_model_ver = 'maxH_default';
    end
    segnames_fold = fullfile(mat_dir, 'cell_data', nrn_model_ver, 'secnames');
    segnames_file = fullfile(segnames_fold, sprintf('secnames%g.mat', cell_id));
    segnames = load(segnames_file);
    segname = segnames.secnames;
    subcellular_names = ["soma" "axon" "Node" "Myelin" "Unmyelin" "dend" "apic"]';
    segindex_soma = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(1)));
    segindex_axon = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(2)));
    segindex_axon = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(3))) + segindex_axon;
    segindex_axon = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(4))) + segindex_axon;
    segindex_axon = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(5))) + segindex_axon;
    segindex_axon = logical(segindex_axon);
    segindex_basal = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(6)));
    segindex_apic = cellfun(@(x) ~isempty(x),strfind(segname,subcellular_names(7)));

end