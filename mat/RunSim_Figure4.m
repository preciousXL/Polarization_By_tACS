mat_dir = addPaths();
nrn_dir = [mat_dir '/../nrn'];
cell_model_names = cellModelNames;
list_cell_id = 1:1:25;
list_f = 10;

polarize_amp = cell(length(list_cell_id), 1);
polarize_phase = cell(length(list_cell_id), 1);

for m=1:length(list_cell_id)
    cell_id = list_cell_id(m);
    Edirect = find_maxpolar_Edirect(cell_id, 'soma');
    nrn_model_ver = 'maxH';
    snowp_test.cell_id = cell_id;
    snowp_test.cell_model_name = cell_model_names{cell_id};
    snowp_test = outputMorphParams(nrn_model_ver, snowp_test);
    model_prefix = 'cell_dat';
    model_name = sprintf('%s_%s',model_prefix,nrn_model_ver);
    [simulation_test,cell_model_test,waveform_test,Efield_test] = loadParams(model_prefix,model_name,snowp_test);
    simulation_test.run_name = [snowp_test.cell_model_name '_' model_name];
    Efield_test.theta = Edirect(1, 1);
    Efield_test.phi = Edirect(2, 1);
    run_params_fold = fullfile(nrn_dir,'params',simulation_test.run_name);        
    if exist(run_params_fold,'dir') == 0                                                  
        mkdir(run_params_fold);                                                           
    end
    tmp_data_folder = fullfile(nrn_dir,'tmp',simulation_test.run_name); 
    if exist(tmp_data_folder,'dir') == 0                                                    
        mkdir(tmp_data_folder);                                                    
    end
    
    for i=1:length(list_f)
        simulation = simulation_test;
        cell_model = cell_model_test;
        waveform = waveform_test;
        Efield = Efield_test;
        tACS = struct();
        tACS.dt = 0.025; 
        tACS.amp = 1;
        tACS.freq = list_f(i);

        tACS.dur = round(1000 * 20/tACS.freq); 

        dt = 0.025;
        DEL = 0;
        DUR = tACS.dur;
        tstop = DEL + DUR + 0;
        [tvec, Evec] = tACSwave(dt, DEL, DUR, tstop, tACS, 0);
        writeVectorBin_snowp(run_params_fold, i, tvec, Evec);

        simulation.dt = dt;
        simulation.tstop = tstop;
        simulation.vm_record_mode = 2;  % section center voltage
        simulation.spike_record_mode = 0;
        simulation.numthetaphi = num2str(i);
        filename = sprintf('defineParams%g.hoc', i);
        personal_name = sprintf('amp%.2fmV_f%gHz_theta%g_phi%g_', tACS.amp, tACS.freq, Efield.theta, Efield.phi);
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
        % compute polarization
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
        snowp.numthetaphi = i;
        snowp.iontype = 'Default';
        [polarize, polarphase] = calc_tACS_polarization_phase(snowp);
    end
    polarize_amp{m, 1} = polarize;
    polarize_phase{m, 1} = polarphase;
end
% 数据保存
savefold = fullfile(mat_dir, 'snowp_data');
if ~exist(savefold, "dir")
    mkdir(savefold);
end
savefile = fullfile(savefold, sprintf('polarization_1mV_10Hz_25cells_Default'));
save(savefile, "polarize_amp", "polarize_phase");