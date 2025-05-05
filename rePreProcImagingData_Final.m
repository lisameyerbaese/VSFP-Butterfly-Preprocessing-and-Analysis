%% this function loads in the raw imaging 
% data and runs the preprocessing pipeline on the data again
function rePreProcImagingData_Final
    % load in the excel sheet that gets loaded in for data analysis
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end
    
    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
    
    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % for now delete VSFP 25
    
    saveFilePath = 'X:\labs\keilholz-lab\Lisa\VSFP ButterFly\Data\VSFP_50Hz\';
    
    % loop through each animal to preprocess the needed data
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        seq_info.temp_reg = 1; 
    
        addpath([startFile, '\VSFP ButterFly\RhettData\Shared_voltage_imaging_Code\reg_files'])
        for k = 1 : length(all_trials)
            %get corresponding pupil data
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;
            fileName = [mouse,'_',num2str(date),'_',num2str(trial)];
            out = preProcVSFP10L_Final('VSFP_01A0', date, trial, mouse, seq_info);
            save([saveFilePath,fileName,'.mat'],'out','-v7.3')
        end
    end

end