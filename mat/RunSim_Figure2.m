%% Generating data for Figure2
%% add paths
mat_dir = addPaths();
nrn_dir = [mat_dir '/../nrn'];

%% Define electric filed direction: totoling 266 directions
list_cell_id = 1:1:25;
list_phi = 0:15:360-1;
list_theta = 0:15:180;
list_params = [list_theta(1) 0];
for i=2:length(list_theta)-1
    for j=1:length(list_phi)
        list_params = [list_params; [list_theta(i) list_phi(j)]];
    end
end
list_params = [list_params; [list_theta(end) 0]];

%% Setting parallel pool
numCPUs = feature('numCores');                                                              
fprintf('Number of CPUs available = %g\n',numCPUs);  
if numCPUs > 1 
    if numCPUs > 4
       pc_storage_dir = fullfile(mat_dir,'pc_storage',getenv('SLURM_JOB_ID'));              
       if ~exist(pc_storage_dir, "dir")
           mkdir(pc_storage_dir);
       end
       pc = parcluster('local');
       pc.JobStorageLocation =  pc_storage_dir;
       pc.NumWorkers = numCPUs;
    else
        pc = parcluster('local');
    end
   poolobj = parpool(pc,numCPUs);
end

%% Run simulation: field intensity is 1mV/mm, field frequency is 10Hz.
for k=1:25
    % parameter setting
    cell_model_names = cellModelNames; % cell names
    num_cells = length(cell_model_names);
    nrn_model_ver = 'maxH';
    cell_id = list_cell_id(k);
    snowp_test.cell_id = cell_id;
    snowp_test.cell_model_name = cellModelNames(cell_id);
    snowp_test = outputMorphParams(nrn_model_ver, snowp_test); % cell morphology scalling
    model_prefix = 'cell_dat';
    model_name = sprintf('%s_%s',model_prefix,nrn_model_ver);
    [simulation_test,cell_model_test,waveform,Efield_test] = loadParams(model_prefix,model_name,snowp_test);
    simulation_test.run_name = [snowp_test.cell_model_name '_' model_name];

    % define tACS waveform
    tACS.dt   = 0.025; % ms
    tACS.amp  = 1;     % mV/mm
    tACS.freq = 10;    % Hz
    tACS.dur = round(1000 * 20 / tACS.freq); % ms, no less than 20 cycles

    dt  = 0.025;
    DEL = 0;           % for that Vm has been initiated
    DUR = tACS.dur;
    tstop = DEL + DUR + 0;
    [tvec, Evec] = tACSwave(dt, DEL, DUR, tstop, tACS, 0);
    run_params_fold = fullfile(nrn_dir,'params',simulation_test.run_name);  % parameter fold                       
    if exist(run_params_fold,'dir') == 0                                                    
        mkdir(run_params_fold);                                                           
    end
    tmp_data_folder = fullfile(nrn_dir,'tmp',simulation_test.run_name);     % temp data fold
    if exist(tmp_data_folder,'dir') == 0                                                   
        mkdir(tmp_data_folder);                                                           
    end

    % change simulation parameters
    simulation_test.dt = dt;
    simulation_test.tstop = tstop;
    simulation_test.vm_record_mode = 2;
    simulation_test.spike_record_mode = 0;
    cell_model_test.replace_axon = 0;
    cell_model_test.myelinate_axon = 1;

    % define data vector
    polarize266 = cell(length(list_params), 1);
    polarize_phase = cell(length(list_params), 1);
    tic;
    % run simulation
    parfor m=1:length(list_params)
        writeVectorBin_snowp(run_params_fold, m, tvec, Evec);   % field data
        Efield = Efield_test;
        simulation = simulation_test;
        cell_model = cell_model_test;
        Efield.theta = list_params(m, 1); 
        Efield.phi = list_params(m, 2); 
        simulation.numthetaphi = num2str(m);
        filename = sprintf('defineParams%g.hoc', m);
        personal_name = sprintf('amp%.2fmV_f%gHz_theta%g_phi%g_', tACS.amp, tACS.freq, Efield.theta, Efield.phi);
        simulation.personal_name = personal_name;
        writeParamsHoc_snowp(run_params_fold,simulation,cell_model,Efield,waveform,filename) % simulation data
        
        % simulation using NEURON
        cd(nrn_dir);   
        nrn_file = fullfile(nrn_dir,'init_field_stim.hoc'); 
        commands = sprintf(' -c "run_name=\\"%s"\\" -c "personal_name=\\"%s"\\" -c "numthetaphi=\\"%s"\\" -c "loadFiles()" -c "run_stim()" -c "quit()"',simulation.run_name, simulation.personal_name, simulation.numthetaphi);                 
        options = sprintf('-nopython -NSTACK %g -NFRAME %g ',100000,20000);
        if ispc
            system(['D:\nrn\bin\nrniv.exe -nobanner ' options nrn_file commands]);            
        elseif isunix       
            system(['./special ' options nrn_file commands]);                                  
        end
        cd(mat_dir);

        % compute polarization length and polarization phase
        snowp = struct();
        snowp.dt = dt;
        snowp.DEL = DEL;
        snowp.DUR = DUR;
        snowp.record_dt = simulation.record_dt;
        snowp.personal_name = personal_name;
        snowp.nrn_dir = nrn_dir;
        snowp.mat_dir = mat_dir;
        snowp.deletefile = 1;
        snowp.run_name = simulation.run_name;
        snowp.amp = tACS.amp;
        snowp.freq = tACS.freq;
        snowp.cell_id = cell_id;
        snowp.numthetaphi = m;
        snowp.iontype = 'Default'; % 'Default', 'Exc', or 'Inh'
        [polarize, polarphase] = calc_tACS_polarization_phase(snowp); 
        polarize266{m, 1} = polarize;
        polarize_phase{m, 1} = polarphase;
    end
    toc;
    % save data for Figure 2
    savefold = fullfile(mat_dir, 'snowp_data', sprintf('%gmV_%gHz', tACS.amp, tACS.freq)); % '..\snowp_data\1mV_10Hz'
    if ~exist(savefold, "dir")
        mkdir(savefold);
    end
    savefile = fullfile(savefold, sprintf('polarization_cell%g', cell_id)); % '..\snowp_data\1mV_10Hz\polarization_cell25.mat'
    save(savefile, "polarize266","polarize_phase");
end