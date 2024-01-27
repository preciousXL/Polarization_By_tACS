%% Setup - run this function once to set up data files necessary for simulations
function InitializeCellData_snowp()
mat_dir = addPaths;   

%% Generate cell data (based on neuron models placed in '../nrn/cells/')
iontype = "Default";
writeCellModelTable(iontype);                                               % creates cell_data/cell_data.mat, called by cellModelNames, "Default", "Exc", or "Inh", same for each iontype
cell_model_names = cellModelNames;                                          % 25*1 cell, outputs all cell_model_names 
num_cells = length(cell_model_names);                                        
nrn_model_ver = 'maxH';                                                     % adult, human, myelinated (see outputMorphParams.m)
for cell_id = 1:num_cells
    coordinates_addr = [mat_dir filesep 'cell_data' filesep sprintf('maxH_%s', iontype) filesep 'coordinates' filesep sprintf('coordinates%g.mat', cell_id)];
    if ~exist(coordinates_addr, "file")
        saveCellData(cell_id,nrn_model_ver,iontype);                                % save coordinates, diameters, areas, secnames ... for each cell
    end
end                    
end