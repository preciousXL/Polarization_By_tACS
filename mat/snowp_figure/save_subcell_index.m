%% find the index of soma, basal, apic, and axon
subcell_index = cell(25, 4);
mat_dir = addPaths;
for m=1:25
    cell_id = m;
    [idxsoma, idxaxon, idxbasal, idxapic] = find_subcell_index(cell_id, mat_dir, 'maxH_default');
    tempindex = (1:1:length(idxsoma))';
    indexsoma = tempindex(idxsoma) - 1; % python index: start with 0
    indexaxon = tempindex(idxaxon) - 1;
    indexbasal = tempindex(idxbasal) - 1;
    indexapic = tempindex(idxapic) - 1;
    subcell_index{m, 1} = indexsoma;
    subcell_index{m, 2} = indexaxon;
    subcell_index{m, 3} = indexbasal;
    subcell_index{m, 4} = indexapic;
end
savefile = sprintf('%s', mat_dir, filesep, 'snowp_data', filesep, '1mV_10Hz_subcell_index');
save(savefile, "subcell_index")