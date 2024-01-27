function Sec_restVm = calc_section_restVm(iontype, list_cell_id)
    if nargin < 2
        list_cell_id = 1:1:25;
    end
    mat_dir = addPaths;
    nrn_dir = [mat_dir '/../nrn'];
    cell_model_names = cellModelNames;

    Sec_restVm = cell(length(list_cell_id), 1);
    parfor i=1:length(list_cell_id)
        cell_id = list_cell_id(i);
        nrn_model_ver = 'maxH';
        snowp1 = struct();
        snowp1.cell_id = cell_id;
        snowp1.cell_model_name = cell_model_names{cell_id};
        snowp1 = outputMorphParams(nrn_model_ver, snowp1);
        model_prefix = 'cell_dat';
        model_name = sprintf('%s_%s',model_prefix,nrn_model_ver);
        [simulation,cell_model,waveform,Efield] = loadParams(model_prefix,model_name,snowp1);
        simulation.run_name = [snowp1.cell_model_name '_' model_name];
        % tACS waveform
        tACS = struct();
        tACS.dt = 0.025; % ms
        tACS.amp = 0; % mV/mm
        tACS.freq = 1; % Hz
        tACS.dur = 0; % ms
        
        dt = 0.025;
        DEL = 10; 
        DUR = tACS.dur;
        tstop = DEL + DUR;
        [tvec, Evec] = tACSwave(dt, DEL, DUR, tstop, tACS, 0);
        run_params_fold = fullfile(nrn_dir,'params',simulation.run_name);                       
        if exist(run_params_fold,'dir') == 0                                                    
            mkdir(run_params_fold);                                                             % './nrn/params/L1_NGC-DA_bNAC219_1_cell_dat_maxH'
        end
        tmp_data_folder = fullfile(nrn_dir,'tmp',simulation.run_name);
        if exist(tmp_data_folder,'dir') == 0                                                    
            mkdir(tmp_data_folder);                                                             % './nrn/params/L1_NGC-DA_bNAC219_1_cell_dat_maxH'
        end
        writeVectorBin_snowp(run_params_fold, 0, tvec, Evec);

        simulation.dt = dt;
        simulation.tstop = tstop;
        simulation.vm_record_mode = 2; % record all resting membrane voltages at x=0.5
        simulation.spike_record_mode = 0;
        cell_model.replace_axon = 0;
        cell_model.myelinate_axon = 1;
        
        Efield.theta = 0;
        Efield.phi = 0;
        simulation.numthetaphi = '0';
        filename = sprintf('defineParams%s.hoc', simulation.numthetaphi);
        personal_name = 'resting_potential_';
        simulation.personal_name = personal_name;
        writeParamsHoc_snowp(run_params_fold,simulation,cell_model,Efield,waveform,filename)
        
        % simulation
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

        % extract resting potential of each sections of each neurons
        file_tvec = fullfile(nrn_dir, 'tmp', simulation.run_name, [personal_name 'tvec.txt']);
        file_restVm = fullfile(nrn_dir, 'tmp', simulation.run_name, [personal_name 'vm_data_bin.txt']);
        fid = fopen(file_restVm, 'r');
        fread(fid, 1, 'double'); 
        lentvec = fread(fid, 1, 'double');
        numSect = fread(fid, 1, 'double');
        vm = zeros(lentvec, numSect);
        for j=1:numSect
            if j ~= 1
                fread(fid, 1, 'double'); 
            end
            vm(:, j) = fread(fid, lentvec, 'double'); 
        end
        fclose('all');
        Vm_rest = vm(end, :)';
        Sec_restVm{i, 1} = Vm_rest;
        % delete files
        delete(file_restVm, file_tvec); 
    end
    % save resting potential file
    savefile = fullfile(mat_dir, 'snowp_data', sprintf('Sec_restVm_%s_25cells.mat', iontype));
    save(savefile, "Sec_restVm");
end


