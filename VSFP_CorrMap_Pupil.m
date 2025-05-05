% Function takes xlsx sheet and generates figures for following analysis
%   1) take pupil diameter and obtains correlation map for volt and hemo
%      (low pass <1Hz) signals
%
% Inputs: 
%        optional mice number, or list of numbers i.e.
%        VSFP_CorrMap_Pupil('VSFP24','VSFP25') or VSFP_CorrMap_Pupil('VSFP24')
%        if no input is given, function runs through all animals
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP_CorrMap_Pupil(varargin)

    
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % delete VSFP 25, no pupil data for these sessions!
    
    if isempty(varargin) == 0
        all_mice = varargin;
    end
    
    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil\'];   % Your destination folder

    % select image used for mask
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    
    % load in cortical map Cortical_Map.png
    disp('Selecting Cortical Map to Use for Registration')
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));
    
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        avgTrial_Hemo = zeros(100,100,length(all_trials));
        avgTrial_Volt = zeros(100,100,length(all_trials));
        
        
        for k = 1:length(all_trials)
            %get corresponding pupil data
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;

            % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            %image_data = strcat('X:\keilholz-lab\Lisa\VSFP ButterFly\RhettData\vsfp_50Hz_proc\',image_file,'.mat');
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            %% Get the Correlation maps 
            %taken from Rhetts code CorrelationMap.m to register the image to Allen

            % get pupil data
            pupil_data2x = vsfp_data.pupil_data2x; 
            
            %voltage signal 
            data_temp = vsfp_data.out.imgDR3;

            % Reshape data temp to 2d
            sXY = vsfp_data.out.sX * vsfp_data.out.sY;
            data_temp_2d = reshape(data_temp,[sXY,vsfp_data.out.sZ]);

            % Correlation With Pupil At Baseline
            corr_map_temp = fastCorr2(data_temp_2d(:,1:length(pupil_data2x))',pupil_data2x');
            corr_map = reshape(corr_map_temp,[100,100]);
            corr_map_reg_volt = corr_map;

            % hemo signal 
            data_temp = -1 * vsfp_data.out.hemoLP;
            
            % Reshape data temp to 2d
            sXY = vsfp_data.out.sX * vsfp_data.out.sY;
            data_temp_2d = reshape(data_temp,[sXY,vsfp_data.out.sZ]);

            % lowpass filter hemo signal at 1 Hz, lowpass fitlers each
            % column independently 
            % fs = vsfp_data.out.Fs;             
            % data_temp = lowpass(data_temp_2d',1,fs,ImpulseResponse="iir",Steepness=0.95);
            % data_temp_2d = data_temp';
            % 
            % Correlation With Pupil At Baseline
            corr_map_temp = fastCorr2(data_temp_2d(:,1:length(pupil_data2x))',pupil_data2x');
            corr_map = reshape(corr_map_temp,[100,100]);
            corr_map_reg_hemo = corr_map;
            
            %% mask out non-cortex
            % get correct pixel size
            maskImg = imresize(maskImg, [vsfp_data.out.sX, vsfp_data.out.sY]); 
            maskImg = maskImg(:,:,1);

            %apply mask on transformed data 
            masked_corr_map_reg_hemo = imrotate(corr_map_reg_hemo, 270);
            maskk = maskImg < 100;
            masked_corr_map_reg_hemo(maskk)= 0;
            
            % get atlas to overlay
            fixed_map = imresize(fixed_map, [vsfp_data.out.sX, vsfp_data.out.sY]);
            M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts

            % plot correlation map for hemo signal
            f1 = figure(1); imagesc(masked_corr_map_reg_hemo), colorbar
            hold on
            a = imagesc(M);
            hold off
            alpha(a,.2)
            title('Max Corr Coef Hemo and Pupil')
            avgTrial_Hemo(:,:,k) = masked_corr_map_reg_hemo;
            
            %apply mask on transformed data 
            masked_corr_map_reg_volt = imrotate(corr_map_reg_volt, 270);
            maskk = maskImg < 100;
            masked_corr_map_reg_volt(maskk)= 0;
            
            % plot correlation map for volt signal
            f2 = figure(2); imagesc(masked_corr_map_reg_volt), colorbar
            hold on
            a = imagesc(M);
            hold off
            alpha(a,.2)
            title('Max Corr Coef Volt and Pupil')
            axis off;
            avgTrial_Volt(:,:,k) = masked_corr_map_reg_volt;
            
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_CorrMap_Hemo'];
            saveFig(f1, fname, FolderName);
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_CorrMap_Volt'];
            saveFig(f2, fname, FolderName);
            close (f1,f2);
        end
        
        avg_Hemo = mean(avgTrial_Hemo,3);
        avg_Volt = mean(avgTrial_Volt,3);
        
        % plot correlation map for avg hemo signal
        f1 = figure(1); imagesc(avg_Hemo), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Max Corr Coef Hemo and Pupil')
        
        % plot correlation map for volt signal
        f2 = figure(2); imagesc(avg_Volt), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Max Corr Coef Volt and Pupil')
        axis off;

        fname = [mouse,'_','avg','_CorrMap_Volt.fig'];
        saveFig(f2, fname, FolderName);
        fname = [mouse,'_','avg','_CorrMap_Hemo.fig'];
        saveFig(f1, fname, FolderName);
        fname = [mouse,'_','avgCorrMap.mat'];
        save(fullfile(FolderName, fname), 'avgTrial_Volt', 'avgTrial_Hemo'); 
        close all;
        
    end   
end