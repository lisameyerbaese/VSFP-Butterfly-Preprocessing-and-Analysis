% Function takes xlsx sheet and generates figures for following analysis
%   1) Randomly selects pupil and imaging data 100 times for shuffled
%   control 
%
% Inputs: 
%        none
%
% Outputs: 
%        .mat files that contains corr values for 100 runs
%
% Written by Lisa Meyer-Baese

function VSFP_CorrMap_Pupil_shuffled(varargin)

   close all
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');

    folderName = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil_shuffled\'];   % Your destination folder
    if ~exist(folderName, 'dir')
       mkdir(folderName)
    end
    
    % select image used for mask
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    
    % load in cortical map Cortical_Map.png
    disp('Selecting Cortical Map to Use for Registration')
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));
    
    % Get list of all good trials to pull pupil data from 
    all_mice = {'VSFP24', 'VSFP27', 'VSFP28', 'VSFP29', 'VSFP30'};
    all_trials_shuffle = [];
    all_dates_shuffle = [];
    for s = 1:length(all_mice)
        temp = find(contains(T.Mouse,all_mice{s}));
        all_trials_shuffle = [all_trials_shuffle;T.Trials(temp)];
        all_dates_shuffle = [all_dates_shuffle; T.Date(temp)];
    end
    
    % loop through and randomly select imaging and pupil dataset 100 times 
    for ind = 1:100
        % randomly pick a trial from the list of all good trials load the
        % pupil trace from this trial 
        rand_trial = randsample( length(all_trials_shuffle),1);
        %get corresponding pupil data for new trial and date
        FindTable = T((T.Trials == all_trials_shuffle(rand_trial) & (T.Date == all_dates_shuffle(rand_trial))),:);
        date = FindTable.Date;
        trial = FindTable.Trials;
        mouse = FindTable.Mouse; 
        mouse = mouse{1};
        
        % load the .mat file that contains the data
        image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
        image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
        vsfp_data=load(image_data);

        % get pupil data
        pupil_data2x = zscore(vsfp_data.pupil_data2x); 
        clear vsfp_data
        
        %% imaging data 
        % randomly pick another trial and load the voltage and hemodynamics imaging data from
        % it
        rand_trial = randsample( length(all_trials_shuffle),1);
        %get corresponding volt/hemo data for new trial and date
        FindTable = T((T.Trials == all_trials_shuffle(rand_trial) & (T.Date == all_dates_shuffle(rand_trial))),:);
        date = FindTable.Date;
        trial = FindTable.Trials;
        mouse = FindTable.Mouse; 
        mouse = mouse{1};
        
        % load the imaging data
        image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
        image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
        vsfp_data=load(image_data);
        
        % Volt and LP Hemo, rotate images to match orientation the allen atlas mask 
        data_temp_volt = vsfp_data.out.imgDR3;
        data_temp_hemo = -1 * vsfp_data.out.hemoLP;
        
        %account for differences in length 
        if length(pupil_data2x) > length(data_temp_volt)
                pupil_data2x = pupil_data2x(1:length(data_temp_volt));
                len_Z = length(data_temp_volt);
        else 
            data_temp_volt = data_temp_volt(:,:,1:length(pupil_data2x));
            data_temp_hemo = data_temp_hemo(:,:,1:length(pupil_data2x));
            len_Z = length(pupil_data2x);
        end
        
        % correlation with voltage data
        % Reshape data temp to 2d
        sXY = vsfp_data.out.sX * vsfp_data.out.sY;
        data_temp_2d = reshape(data_temp_volt,[sXY,len_Z ]);

        % Correlation With Pupil At Baseline
        corr_map_temp = fastCorr2(data_temp_2d(:,1:length(pupil_data2x))',pupil_data2x');
        corr_map_reg_volt = reshape(corr_map_temp,[100,100]);

        % mask out non-cortex
        maskImg = imresize(maskImg, [vsfp_data.out.sX, vsfp_data.out.sY]); 
        maskImg = maskImg(:,:,1);
        
        % % Register to mouse template
        % image_data = strcat('X:\keilholz-lab\Lisa\VSFP ButterFly\RhettData\Shared_voltage_imaging_Code\reg_files\',mouse);
        % temp_load = load([image_data '_reg.mat']);
        % temp_base = temp_load.reg.baseImg;
        % [optimizer, metric] = imregconfig('monomodal'); % monomodal when all cameras are donor channel
        % tform_tempReg = imregtform(vsfp_data.out.baseD,temp_base,'rigid',optimizer,metric);
        % map_ref = imref2d([100,100]);
        % corr_map_reg_volt = imwarp(corr_map(:,:),tform_tempReg,'nearest','OutputView',map_ref);
        
        
        %apply mask on transformed data 
        masked_corr_map_reg_volt = imrotate(corr_map_reg_volt, 270);
        maskk = maskImg < 100;
        masked_corr_map_reg_volt(maskk)= 0;
        
        % hemo signal 
        sXY = vsfp_data.out.sX * vsfp_data.out.sY;
        data_temp_2d = reshape(data_temp_hemo,[sXY,len_Z ]);
        
        % low pass hemo to 1Hz 
        fs = vsfp_data.out.Fs;             
        data_temp = lowpass(data_temp_2d',1,fs,ImpulseResponse="iir",Steepness=0.95);
        data_temp_2d = data_temp';
        
        % Correlation With Pupil At Baseline
        corr_map_temp = fastCorr2(data_temp_2d(:,1:length(pupil_data2x))',pupil_data2x');
        corr_map_reg_hemo = reshape(corr_map_temp,[100,100]);
       
        % % Register to mouse template
        % image_data = strcat('X:\keilholz-lab\Lisa\VSFP ButterFly\RhettData\Shared_voltage_imaging_Code\reg_files\',mouse);
        % temp_load = load([image_data '_reg.mat']);
        % temp_base = temp_load.reg.baseImg;
        % [optimizer, metric] = imregconfig('monomodal'); % monomodal when all cameras are donor channel
        % tform_tempReg = imregtform(vsfp_data.out.baseD,temp_base,'rigid',optimizer,metric);
        % map_ref = imref2d([100,100]);
        % corr_map_reg_hemo = imwarp(corr_map(:,:),tform_tempReg,'nearest','OutputView',map_ref);

        %apply mask on transformed data 
        masked_corr_map_reg_hemo = imrotate(corr_map_reg_hemo, 270);
        maskk = maskImg < 100;
        masked_corr_map_reg_hemo(maskk)= 0;

        % get atlas to overlay
        fixed_map = imresize(fixed_map, [vsfp_data.out.sX, vsfp_data.out.sY]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
        
        avgTrial_Volt(:,:,ind) = masked_corr_map_reg_volt;
        avgTrial_Hemo(:,:,ind) = masked_corr_map_reg_hemo;
        
    end
    
    save([folderName,'AvgTrial_Volt.mat'], 'avgTrial_Volt')
    save([folderName,'AvgTrial_Hemo.mat'], 'avgTrial_Hemo')
end 
    
