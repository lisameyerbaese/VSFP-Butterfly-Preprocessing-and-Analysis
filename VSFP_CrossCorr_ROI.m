% Function takes xlsx sheet and generates figures for following analysis
%   New version of VSFP_Lag_ROI (just a more intuitive function name)
%   calculates the lag between voltage and hemo, hemo and pupil, and volt and pupil per selected ROIS 
%
%
% Inputs: 
%        regionNum = ROI number to select 1-4;
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP_CrossCorr_ROI(regionNum)
    
    close all;
    
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    % file saving destination
    FolderName = [startFile,'\VSFP ButterFly\Data\VSFP_CrossCorr_ROI\', num2str(regionNum)]; 
    %check to see if Folder Exists to Save data, if not create it
    if ~exist(FolderName, 'dir')
      mkdir(FolderName)
    end

    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % delete VSFP 25, no pupil data for these sessions!
    
    % % load in cortical map Cortical_Map.png
    % disp('Selecting Cortical Map to Use for Registration')
    % file = 'Cortical_Map.png';
    % disp('Selecting Cortical Map to Use for Masking')
    % path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    % [fixed_map,~] = imread(fullfile(path_map, file));
    % fixed_map = imresize(fixed_map, [100, 100]);
    % %M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
    
    for roi = regionNum:regionNum
        
        lagAll_volt_pupil = [];
        lagAll_hemo_pupil = [];
        corrAll_volt_pupil = [];
        corrAll_hemo_pupil = [];
        corrAll_volt_face = [];
        corrAll_hemo_face = [];
        lagAll_volt_face = [];
        lagAll_hemo_face = [];
        
        ROI = num2str(roi);
        for m = 1: length(all_mice) 

            mouse = all_mice{m};
            ind = find(contains(T.Mouse,mouse));
            all_trials = T.Trials(ind);
            all_dates = T.Date(ind);

            disp('Checking to see if .mat Avg Corr Map already exists')
            fname = [mouse,'_','avgCorrMap.mat'];
            dataFolder  = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil'];
            if isfile(fullfile(dataFolder, fname))
               disp('File Exists!')
               load(fullfile(dataFolder, fname), 'avgTrial_Volt');
            else 
               disp('File does not exists, running VSFP_CorrMap_Pupil.m to generate needed file')
               VSFP_CorrMap_Pupil(mouse)
               load(fullfile(dataFolder, fname), 'avgTrial_Volt');
            end

            f1 = figure(1);
            disp('Select desired ROI from cortical map')
            avgTrial_Volt = mean(avgTrial_Volt,3);
            imagesc(avgTrial_Volt), title('Select ROI: Avg Corr Map for Volt and Pupil'), colorbar
            r1 = drawcircle('Label','ROI','Color',[1 0 0]);
            maskImg = createMask(r1);

          
            lag_all_volt_hemo = [];
            corr_all_volt_hemo = [];
            lag_all_hemo_pupil = [];
            corr_all_hemo_pupil = [];
            lag_all_volt_pupil = [];
            corr_all_volt_pupil = [];
            lag_all_hemo_face = [];
            corr_all_hemo_face = [];
            lag_all_volt_face = [];
            corr_all_volt_face = [];

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
               
                % get pupil data
                pupil_data2x = zscore(vsfp_data.pupil_data2x); 

                % get orofacial data too 
                face_data2x = vsfp_data.face_data2x'; 
                
                %voltage signal 
                data_volt = vsfp_data.out.imgDR3;

                %hemo signal
                data_hemo = -1 * vsfp_data.out.hemoLP;

                %apply mask on volt and hemo data 
                masked_volt = imrotate(data_volt, 270) .* maskImg;
                masked_hemo = imrotate(data_hemo, 270).* maskImg;

                % Find the regional "Global Signal"
                volt2D = reshape(masked_volt, 10000,[]);
                voltTemp = zscore(volt2D,[],2);
                gsVolt = mean(voltTemp,1);

                hemo2D = reshape(masked_hemo, 10000,[]);
                hemoTemp = zscore(hemo2D,[],2);

                Fs = vsfp_data.out.Fs;  % Sampling Frequency (Hz)

                % get the "global" hemo signal and low pass filter it
                gsHemo = mean(hemoTemp,1);
                gsHemo = lowpass(gsHemo,1,Fs,ImpulseResponse="iir",Steepness=0.95);

                % Cross Corr hemo and volt
                [r,lags] = xcorr(gsHemo, gsVolt, 500,'normalized');
                t = lags/Fs;
                len = round(length(t)/ 2);
                ind = (len - 250: len + 250);
                x = t(ind);
                y = r(ind);
                
                f2 = figure(2);
                subplot(3,1,1)
                plot(gsHemo,'DisplayName','Avg Global Hemo')
                hold on
                plot(gsVolt, 'DisplayName','Avg Global Volt')
                legend
                subplot(3,1,2)
                plot(pupil_data2x,'DisplayName','Pupil Trace')
                subplot(3,1,3)
                plot(face_data2x,'DisplayName','Face Trace')
              
               
                % save figure
                fname = [mouse,'_',num2str(trial),ROI,'Avg Signal in ROI.fig'];
                saveFig(f2, fname, FolderName)
                close (f2)

                f3 = figure(3);
                plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
                hold on
                title(['Response from ROI',ROI,'Cross Corr Hemo/Volt'])
                lag_all_volt_hemo = [lag_all_volt_hemo; x];
                corr_all_volt_hemo = [corr_all_volt_hemo;y];

                % Cross Corr Volt and Pupil
                [r,lags] = xcorr(gsVolt(1:length(pupil_data2x)), pupil_data2x, 500,'normalized');
                t = lags/Fs;
                len = round(length(t)/ 2);
                ind = (len - 250: len + 250);
                x = t(ind);
                y = r(ind);

                f4 = figure(4);
                plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
                hold on
                title(['Response from ROI ',ROI,' Volt and Pupil'])
                lag_all_volt_pupil = [lag_all_volt_pupil; x];
                corr_all_volt_pupil = [corr_all_volt_pupil;y];

                %Cross Corr Hemo and Pupil
                [r,lags] = xcorr(gsHemo(1:length(pupil_data2x)), pupil_data2x, 500,'normalized');
                t = lags/Fs;
                len = round(length(t)/ 2);
                ind = (len - 250: len + 250);
                x = t(ind);
                y = r(ind);

                f5 = figure(5);
                plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
                hold on
                title(['Response from ROI',ROI,'Hemo and Pupil'])
                lag_all_hemo_pupil = [lag_all_hemo_pupil; x];
                corr_all_hemo_pupil = [corr_all_hemo_pupil;y];

                % Cross Corr Volt and Face
                [r,lags] = xcorr(gsVolt(1:length(face_data2x)), face_data2x, 500,'normalized');
                t = lags/Fs;
                len = round(length(t)/ 2);
                ind = (len - 250: len + 250);
                x = t(ind);
                y = r(ind);

                f6 = figure(6);
                plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
                hold on
                title(['Response from ROI ',ROI,' Volt and Face'])
                lag_all_volt_face = [lag_all_volt_face; x];
                corr_all_volt_face = [corr_all_volt_face;y];

                %Cross Corr Hemo and Face
                [r,lags] = xcorr(gsHemo(1:length(face_data2x)), face_data2x, 500,'normalized');
                t = lags/Fs;
                len = round(length(t)/ 2);
                ind = (len - 250: len + 250);
                x = t(ind);
                y = r(ind);

                f7 = figure(7);
                plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
                hold on
                title(['Response from ROI',ROI,'Hemo and Face'])
                lag_all_hemo_face = [lag_all_hemo_face; x];
                corr_all_hemo_face = [corr_all_hemo_face;y];

            end

            fname = [mouse,'ROI',ROI,'.fig'];
            saveFig(f1, fname, FolderName);
            fname = [mouse,ROI,'Corr Volt and Hemo.fig'];
            saveFig(f3, fname, FolderName);

            % calc standard error and plot that for each animal
            err = std(corr_all_volt_hemo,[],1)/sqrt(size(corr_all_volt_hemo,1));
            f8 = figure(8);
            shade(mean(lag_all_volt_hemo,1), mean(corr_all_volt_hemo,1)+err,mean(lag_all_volt_hemo,1),mean(corr_all_volt_hemo,1)- err, 'FillType',[1 2])
            hold on 
            plot(mean(lag_all_volt_hemo,1), mean(corr_all_volt_hemo,1), 'k', 'LineWidth',2)
            title(['Response from ROI ',ROI,' Hemo and Volt'])
            hold off
            
            fname = ['Avg_',mouse,ROI,'Corr Hemo and Volt.fig'];
            saveFig(f8, fname, FolderName);
            fname = [mouse,ROI,'Corr Volt and Pupil.fig'];
            saveFig(f4, fname, FolderName);
            err = std(corr_all_volt_pupil,[],1)/sqrt(size(corr_all_volt_pupil,1));
            
            f9 = figure(9);
            shade(mean(lag_all_volt_pupil,1), mean(corr_all_volt_pupil,1)+err,mean(lag_all_volt_pupil,1),mean(corr_all_volt_pupil,1)- err, 'FillType',[1 2])
            hold on 
            plot(mean(lag_all_volt_pupil,1), mean(corr_all_volt_pupil,1), 'k', 'LineWidth',2)
            title(['Response from ROI ',ROI,' Volt and Pupil'])
            hold off
            fname = ['Avg_',mouse,ROI,'Corr Volt and Pupil.fig'];
            saveFig(f9, fname, FolderName);

            fname = [mouse,ROI,'Corr Hemo and Pupil.fig'];
            saveFig(f5, fname, FolderName);
            err = std(corr_all_hemo_pupil,[],1)/sqrt(size(corr_all_hemo_pupil,1));
            
            f10 = figure(10);
            shade(mean(lag_all_hemo_pupil,1), mean(corr_all_hemo_pupil,1)+err,mean(lag_all_hemo_pupil,1),mean(corr_all_hemo_pupil,1)- err, 'FillType',[1 2])
            hold on 
            plot(mean(lag_all_hemo_pupil,1), mean(corr_all_hemo_pupil,1), 'k', 'LineWidth',2)
            title(['Response from ROI ',ROI,' Hemo and Pupil'])
            hold off
            fname = ['Avg_',mouse,ROI,'Corr Hemo and Pupil.fig'];
            saveFig(f10, fname, FolderName);

            % now for the orofacial movement data
            fname = [mouse,ROI,'Corr Volt and Face.fig'];
            saveFig(f6, fname, FolderName);
            err = std(corr_all_volt_face,[],1)/sqrt(size(corr_all_volt_face,1));
            
            f11 = figure(11);
            shade(mean(lag_all_volt_face,1), mean(corr_all_volt_face,1)+err,mean(lag_all_volt_face,1),mean(corr_all_volt_face,1)- err, 'FillType',[1 2])
            hold on 
            plot(mean(lag_all_volt_face,1), mean(corr_all_volt_face,1), 'k', 'LineWidth',2)
            title(['Response from ROI ',ROI,' Volt and Face'])
            hold off
            fname = ['Avg_',mouse,ROI,'Corr Volt and Face.fig'];
            saveFig(f11, fname, FolderName);

            fname = [mouse,ROI,'Corr Volt and Face.fig'];
            saveFig(f7, fname, FolderName);
            err = std(corr_all_volt_face,[],1)/sqrt(size(corr_all_volt_face,1));
            
            f12 = figure(12);
            shade(mean(lag_all_hemo_face,1), mean(corr_all_hemo_face,1)+err,mean(lag_all_hemo_face,1),mean(corr_all_hemo_face,1)- err, 'FillType',[1 2])
            hold on 
            plot(mean(lag_all_hemo_face,1), mean(corr_all_hemo_face,1), 'k', 'LineWidth',2)
            title(['Response from ROI ',ROI,' Hemo and Face'])
            hold off
            fname = ['Avg_',mouse,ROI,'Corr Hemo and Face.fig'];
            saveFig(f12, fname, FolderName);
            close all;

            lagAll_volt_pupil = [lagAll_volt_pupil; lag_all_volt_pupil ];
            corrAll_volt_pupil = [corrAll_volt_pupil; corr_all_volt_pupil ];
            lagAll_hemo_pupil = [lagAll_hemo_pupil; lag_all_hemo_pupil];
            corrAll_hemo_pupil = [corrAll_hemo_pupil; corr_all_hemo_pupil ];
            lagAll_volt_face = [lagAll_volt_face; lag_all_volt_face ];
            corrAll_volt_face = [corrAll_volt_face; corr_all_volt_face ];
            lagAll_hemo_face = [lagAll_hemo_face; lag_all_hemo_face];
            corrAll_hemo_face = [corrAll_hemo_face; corr_all_hemo_face ];
        end  
        
        % avg for pupil data
        err = std(corrAll_volt_pupil,[],1)/sqrt(size(corrAll_volt_pupil,1));
        f13 = figure(13);
        shade(mean(lagAll_volt_pupil,1), mean(corrAll_volt_pupil,1)+err,mean(lagAll_volt_pupil,1),mean(corrAll_volt_pupil,1)- err, 'FillType',[1 2])
        hold on 
        plot(mean(lagAll_volt_pupil,1), mean(corrAll_volt_pupil,1), 'k', 'LineWidth',2)
        title(['Avg Response from ROI ',ROI,' Volt and Pupil'])
        hold off
        fname = ['Avg_',ROI,'Corr Volt and Pupil.fig'];
        saveFig(f13, fname, FolderName);

        err = std(corrAll_hemo_pupil,[],1)/sqrt(size(corrAll_hemo_pupil,1));
        f14 = figure(14);
        shade(mean(lagAll_hemo_pupil,1), mean(corrAll_hemo_pupil,1)+err,mean(lagAll_hemo_pupil,1),mean(corrAll_hemo_pupil,1)- err, 'FillType',[1 2])
        hold on 
        plot(mean(lagAll_hemo_pupil,1), mean(corrAll_hemo_pupil,1), 'k', 'LineWidth',2)
        title(['Avg Response from ROI ',ROI,' Hemo and Pupil'])
        hold off
        fname = ['Avg_',ROI,'Corr Hemo and Pupil.fig'];
        saveFig(f14, fname, FolderName);
        
        % now for face data
        err = std(corrAll_volt_face,[],1)/sqrt(size(corrAll_volt_face,1));
        f15 = figure(15);
        shade(mean(lagAll_volt_face,1), mean(corrAll_volt_face,1)+err,mean(lagAll_volt_face,1),mean(corrAll_volt_face,1)- err, 'FillType',[1 2])
        hold on 
        plot(mean(lagAll_volt_face,1), mean(corrAll_volt_face,1), 'k', 'LineWidth',2)
        title(['Avg Response from ROI ',ROI,' Volt and Face'])
        hold off
        fname = ['Avg_',ROI,'Corr Volt and Face.fig'];
        saveFig(f15, fname, FolderName);

        err = std(corrAll_hemo_face,[],1)/sqrt(size(corrAll_hemo_face,1));
        f16 = figure(16);
        shade(mean(lagAll_hemo_face,1), mean(corrAll_hemo_face,1)+err,mean(lagAll_hemo_face,1),mean(corrAll_hemo_face,1)- err, 'FillType',[1 2])
        hold on 
        plot(mean(lagAll_hemo_face,1), mean(corrAll_hemo_face,1), 'k', 'LineWidth',2)
        title(['Avg Response from ROI ',ROI,' Hemo and Face'])
        hold off
        fname = ['Avg_',ROI,'Corr Hemo and Face.fig'];
        saveFig(f16, fname, FolderName);
        
    end
end