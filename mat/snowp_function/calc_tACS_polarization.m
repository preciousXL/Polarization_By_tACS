function polarize = calc_tACS_polarization(instruct, methods)
    if nargin < 2
        methods = 2;
    end
    record_dt = instruct.record_dt;
    DEL = instruct.DEL;
    DUR = instruct.DUR;
    personal_name = instruct.personal_name;
    nrn_dir = instruct.nrn_dir;
    mat_dir = instruct.mat_dir;
    deletefile = instruct.deletefile;
    run_name = instruct.run_name;
    amp = instruct.amp;
    cell_id = instruct.cell_id;
    iontype = instruct.iontype;

    Vmrest_file = fullfile(mat_dir, sprintf("snowp_data/Sec_restVm_%s_25cells.mat", iontype));
    if ~exist(Vmrest_file, "file")
        calc_section_restVm(iontype);
    end
    data = load(Vmrest_file);
    Sec_restVm_25cells = data.Sec_restVm;
    sec_restVm = Sec_restVm_25cells{cell_id};
    file_tvec = fullfile(nrn_dir, 'tmp', run_name, [personal_name 'tvec.txt']);
    file_vm = fullfile(nrn_dir, 'tmp', run_name, [personal_name 'vm_data_bin.txt']);
    fid = fopen(file_vm, 'r');
    fread(fid, 1, 'double'); 
    lentvec = fread(fid, 1, 'double');
    numSect = fread(fid, 1, 'double');
    vm = zeros(lentvec, numSect);
    for i=1:numSect
        if i ~= 1
            fread(fid, 1, 'double'); 
        end
        vm(:, i) = fread(fid, lentvec, 'double'); 
    end
    fclose('all');

    if deletefile == 1
        delete(file_vm, file_tvec)
    end

    tstart = round(DEL/record_dt) + 1;
    tend = round((DEL+DUR)/record_dt);
    % tvec_slice = tvec(tstart:tend, 1);
    vm_slice = vm(tstart:tend, :);
    
    polarize = zeros(numSect, 1);
    idx_max = islocalmax(vm_slice, 1);
    idx_min = islocalmin(vm_slice, 1);
    for i=1:numSect
        tempvm = vm_slice(:, i);
        maxvalue = tempvm(idx_max(:, i));
        minvalue = tempvm(idx_min(:, i));
        if length(maxvalue) ~= length(minvalue)
            templen = min(length(maxvalue), length(minvalue));
            maxvalue = maxvalue(end-templen+1:end, 1);
            minvalue = minvalue(end-templen+1:end, 1);
        end
        if methods==1
            polarize(i, 1) = (mean(maxvalue(3:end, 1)) - mean(minvalue(3:end, 1))) / 2;
        elseif methods==2
            polarize(i, 1) = maxvalue(end) - sec_restVm(i, 1);
        end
    end
    polarize = polarize / amp;
end