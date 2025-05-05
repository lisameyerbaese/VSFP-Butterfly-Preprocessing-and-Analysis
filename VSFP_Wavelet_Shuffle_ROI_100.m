% Function takes xlsx sheet and generates figures for following analysis
%   Similar to VSFP_Wavelet_Shuffle_100.m except that this one performs the
%   coherence analysis for each ROI
%
%   1) Randomly selects pupil and imaging data 100 times for shuffled
%   control (randomly selects different traces also for volt and hemo)
%   2) Used to get shuffled controls for the coherence analysis 
%
% Inputs: 
%        which region to look at, i.e any number 1-4
%
% Outputs: 
%        .mat files that contains corr values for 100 runs
%
% Written by Lisa Meyer-Baese

function VSFP_Wavelet_Shuffle_ROI_100(region)

    close all
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    muHV = [];
    muPV =[];
    muPH = [];
    muFV = [];
    muFH =[];
    
    % Get list of all good trials to pull pupil data from 
    all_mice = {'VSFP24', 'VSFP27', 'VSFP28', 'VSFP29', 'VSFP30'};
    all_trials_shuffle = [];
    all_dates_shuffle = [];
    for s = 1:length(all_mice)
        temp = find(contains(T.Mouse,all_mice{s}));
        all_trials_shuffle = [all_trials_shuffle;T.Trials(temp)];
        all_dates_shuffle = [all_dates_shuffle; T.Date(temp)];
    end
    
    for ROI = region:region
        %select same ROI for all 100 shuffled versions
        %using mouse VSFP24
        disp('Checking to see if .mat Avg Corr Map already exists')
        FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil\'];
        fname = [all_mice{1},'_','avgCorrMap.mat'];
        if isfile(fullfile(FolderName, fname))
           disp('File Exists!')
           load(fullfile(FolderName, fname), 'avgTrial_Volt');
        else 
           disp('File does not exists, running VSFP_CorrMap_Pupil.m to generate needed file')
           VSFP_CorrMap_Pupil(all_mice{1})
           load(fullfile(FolderName, fname), 'avgTrial_Volt');
        end
        disp('Select desired ROI from cortical map')
        avgTrial_Volt = mean(avgTrial_Volt,3);
        imagesc(avgTrial_Volt), title(['Select ROI ', num2str(region)]), colorbar
        r1 = drawcircle('Label','ROI','Color',[1 0 0]);
        maskImg = createMask(r1);
    
        % loop through and randomly select imaging and pupil dataset 100 times 
        for ind = 1:100
            % randomly sample 4 trials to load pupil, face, volt and hemo
            % data from
            rand_trial = randsample( length(all_trials_shuffle),4);

            %get corresponding face data for new trial and date
            FindTable = T((T.Trials == all_trials_shuffle(rand_trial(1)) & (T.Date == all_dates_shuffle(rand_trial(1)))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;
            mouse = FindTable.Mouse; 
            mouse = mouse{1};
            
            % load the .mat file that contains the data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data, 'pupil_data2x');
    
            % get pupil data
            pupil_data2x = zscore(vsfp_data.pupil_data2x); 

            %get face data
            FindTable = T((T.Trials == all_trials_shuffle(rand_trial(2)) & (T.Date == all_dates_shuffle(rand_trial(2)))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;
            mouse = FindTable.Mouse; 
            mouse = mouse{1};
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data, 'face_data2x');
            face_data2x = zscore(vsfp_data.face_data2x)'; 

%             minInd = min(length(face_data2x),length(pupil_data2x));
%             face_data2x = face_data2x(1:minInd);
%             pupil_data2x = pupil_data2x(1:minInd);
            clear vsfp_data
            
            %% volt imaging data 
            %get corresponding volt/hemo data for new trial and date
            FindTable = T((T.Trials == all_trials_shuffle(rand_trial(3)) & (T.Date == all_dates_shuffle(rand_trial(3)))),:);
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
            clear vsfp_data;

            %% randomly pick another trial for Hemo data
            %get corresponding volt/hemo data for new trial and date
            FindTable = T((T.Trials == all_trials_shuffle(rand_trial(4)) & (T.Date == all_dates_shuffle(rand_trial(4)))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;
            mouse = FindTable.Mouse; 
            mouse = mouse{1};
            
            % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);
            data_temp_hemo = -1 * vsfp_data.out.hemoLP;
  
            % account for differences in length
            minInd1 = min(length(face_data2x),length(pupil_data2x));
            minInd2 = min(length(data_temp_hemo),length(data_temp_volt));
            minInd = min(minInd1, minInd2);
            data_temp_hemo = data_temp_hemo(:,:,1:minInd);
            data_temp_volt = data_temp_volt(:,:,1:minInd);
            face_data2x = face_data2x(1:minInd);
            pupil_data2x = pupil_data2x(1:minInd);

            
%             %account for differences in length 
%             if length(face_data2x) > length(data_temp_volt)
%                     face_data2x = face_data2x(1:length(data_temp_volt));
%                     pupil_data2x = pupil_data2x(1:length(data_temp_volt));
%                     len_Z = length(data_temp_volt);
%             else 
%                 data_temp_volt = data_temp_volt(:,:,1:length(face_data2x));
%                 data_temp_hemo = data_temp_hemo(:,:,1:length(face_data2x));
%                 len_Z = length(face_data2x);
%             end
            
            %% mask out all other areas
            % get correct pixel size
            maskImg = imresize(maskImg, [vsfp_data.out.sX, vsfp_data.out.sY]); 
            maskImg = maskImg(:,:,1);
    
            %apply mask on volt and hemo data 
            masked_volt = imrotate(data_temp_volt, 270);
            masked_hemo = imrotate(data_temp_hemo, 270);
            
            masked_volt = masked_volt .* maskImg; 
            masked_hemo = masked_hemo .* maskImg; 
            
            % Find the regional "Global Signal"
            volt2D = reshape(masked_volt, 10000,[]);
            voltTemp = zscore(volt2D,[],2);
            gsVolt = mean(voltTemp,1);
    
            hemo2D = reshape(masked_hemo, 10000,[]);
            hemoTemp = zscore(hemo2D,[],2);
            gsHemo = mean(hemoTemp,1);
    
            avg_global_volt =  gsVolt;
            avg_global_hemo =  gsHemo;
    
            Fs = vsfp_data.out.Fs;             
            t = 0:1/Fs:length(pupil_data2x)/Fs-1/Fs;
    
            %% Calc coherence 
            % take data from 20 to 60 s s:s1
            s = find(t == 20);
            s1 = find( round(t, 2) == 60);
            subplot(5,1,1)
            [wcoh,~,F,coi] =wcoherence(avg_global_hemo(1:length(pupil_data2x)),avg_global_volt(1:length(pupil_data2x)),Fs);
            muHV = [muHV,mean(wcoh(1:110,s:s1),2)];
            subplot(5,1,2)
            [wcoh,~,F,coi] = wcoherence(pupil_data2x,avg_global_volt(1:length(pupil_data2x)),Fs);
            muPV = [muPV,mean(wcoh(1:110,s:s1),2)];
            subplot(5,1,3)
            [wcoh,~,F,coi] = wcoherence(pupil_data2x,avg_global_hemo(1:length(pupil_data2x)),Fs);
            muPH = [muPH,mean(wcoh(1:110,s:s1),2)];
    
            subplot(5,1,4)
            [wcoh,~,F,coi] = wcoherence(face_data2x,avg_global_volt(1:length(face_data2x)),Fs);
            muFV = [muFV,mean(wcoh(1:110,s:s1),2)];
            subplot(5,1,5)
            [wcoh,~,F,coi] = wcoherence(face_data2x,avg_global_hemo(1:length(face_data2x)),Fs);
            muFH = [muFH,mean(wcoh(1:110,s:s1),2)];
    
            
        end
    end
    
    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_Wavelet_Shuffle\', num2str(ROI),'.mat'];
    % save variables depending on which ROI was used
    switch ROI
        case 1
            muHV_1 = muHV;
            muPV_1 = muPV;
            muPH_1 = muPH; 
            muFV_1 = muFV;
            muFH_1 = muFH; 
            save(FolderName, 'muHV_1','muPV_1','muPH_1','muFV_1','muFH_1')      
        case 2
            muHV_2 = muHV;
            muPV_2 = muPV;
            muPH_2 = muPH; 
            muFV_2 = muFV;
            muFH_2 = muFH;
            save(FolderName, 'muHV_2','muPV_2','muPH_2','muFV_2','muFH_2') 
        case 3
            muHV_3 = muHV;
            muPV_3 = muPV;
            muPH_3 = muPH; 
            muFV_3 = muFV;
            muFH_3 = muFH; 
            save(FolderName, 'muHV_3','muPV_3','muPH_3','muFV_3','muFH_3') 
        case 4
            muHV_4 = muHV;
            muPV_4 = muPV;
            muPH_4 = muPH; 
            muFV_4 = muFV;
            muFH_4 = muFH; 
            save(FolderName, 'muHV_4','muPV_4','muPH_4','muFV_4','muFH_4') 
    end
end
    
