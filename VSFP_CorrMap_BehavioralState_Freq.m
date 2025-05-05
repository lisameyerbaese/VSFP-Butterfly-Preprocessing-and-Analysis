% Function takes xlsx sheet and generates figures for following analysis
%   Modified version of VSFP_CorrMap_Face_Freq_Avg that has a more
%   intuative name and loads in the data more efficiently
%
%   NEED TO SPECIFIY  min_sec = 2 PARAMETER AS WELL AS PASSBAND CUTOFF
%
%   Instead of combining all periods of mvmt and rest to do corr, corr is
%   found first and then avg is taken 
%   1) Based on orofacial movements, splits data into rest and mvmt periods
%  
%   2) take pupil diameter and obtains correlation map for hemo and volt
%      signals
%
% Inputs: 
%        min_sec = 2;
%        type of filtering, either BP or LP
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP_CorrMap_BehavioralState_Freq(min_sec, filterType)

    close all;
    
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    %check to see if Folder Exists to Save data, if not create it
    FolderName = [startFile,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\', filterType,'_',num2str(min_sec),'\']; 
    if ~exist(FolderName, 'dir')
      mkdir(FolderName)
    end
   
    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % delete VSFP 25, no pupil data
    Fs = 50;

    % specificy the type of filter to be used 
    if isequal(filterType,'LP')
         %Low Pass Filter 
         dataFilt = make_ChebII_filter(2, Fs, 0.15, 0.2, 20);
    elseif isequal(filterType,'BP') 
        dataFilt = make_ChebII_filter(1, Fs,[0.15 1], [0.1 1.5], 20); % bandpass from 0.15 to 1 Hz
    else 
        error('Incorrect or missing filter type')
        return
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
    
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        avgTrial_Hemo_mvmt = zeros(100,100,length(all_trials));
        avgTrial_Hemo_rest = zeros(100,100,length(all_trials));
        avgTrial_Volt_mvmt = zeros(100,100,length(all_trials));
        avgTrial_Volt_rest = zeros(100,100,length(all_trials));
        
        for k = 1:length(all_trials)
             
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;

            % load the data for this trial
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            % voltage and hemo data
            data_volt = vsfp_data.out.imgDR3;
            data_hemo = - 1 .* vsfp_data.out.hemoLP;

            % filter voltage and hemo data
            for x = 1:100
                for y = 1:100
                    %Filter data 
                    data_volt(x,y,:) = filtfilt(dataFilt.numer, dataFilt.denom,squeeze(squeeze(data_volt(x,y,:)))');
                    data_hemo(x,y,:) = filtfilt(dataFilt.numer, dataFilt.denom,squeeze(squeeze(data_hemo(x,y,:)))');
                end
            end 

            % get pupil data
            pupil_data2x = zscore(vsfp_data.pupil_data2x);

            % get face data
            face_data2x = zscore(vsfp_data.face_data2x);
        
            %mvmt_filt_z = zscore(face_data2x);
            % find periods that are above 0 to classify them as
            % movement/rest
            mvmt_thr = 0; %1.5.*std(mvmt_filt_z);
            mvmt_inds = find( face_data2x>mvmt_thr);
            rest_inds = find( face_data2x<mvmt_thr);

            % Find continuous periods of rest and movement
            mvmt_inds01 = zeros(1,length(mvmt_inds));
            mvmt_inds01(mvmt_inds) = 1;
            mvmt_inds01(1) = 0; mvmt_inds01(end) = 0;
            mvmt_pulses = findPulses(mvmt_inds01);

            rest_inds01 = zeros(1,length(rest_inds));
            rest_inds01(rest_inds) = 1;
            rest_inds01(1) = 0; rest_inds01(end) = 0;
            rest_pulses = findPulses(rest_inds01); 
            
            %duration of mvmt and rest periods
            mvmt_periods = mvmt_pulses.ends-mvmt_pulses.starts;
            rest_periods = rest_pulses.ends-rest_pulses.starts;

            %% Histogram of rest vs movement period lengths
            edges = 1:5:500;
            plot_fig = 1;
            if plot_fig == 1
                f1 = figure(1); 
                histogram(mvmt_periods,edges), hold on
                histogram(rest_periods,edges)
                title(' Histogram of Rest vs Movement Period Lengths')
                legend ('Movement', 'Rest')
            end
        
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_Histogram'];
            saveFig(f1, fname, FolderName);
            close(f1)
            %%

            min_length = min_sec.*Fs;
            mvmt_good = find(mvmt_periods>min_length);
            rest_good = find(rest_periods>min_length);

            %% Sample plot
            plot_fig = 1;
            if plot_fig == 1
                f3 = figure(3); 
                times = 1/Fs:1/Fs:length(face_data2x).*1/Fs;
                plot(times(1:end), face_data2x)
                hold on
                plot(times(1:end), pupil_data2x(1:length(times)))
                line([times(1),times(end)],[mvmt_thr,mvmt_thr])
                scatter(times(mvmt_inds),ones(length(mvmt_inds),1).*max(face_data2x))
                title('Time Course with thresholded movements')
            end
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_Traces'];
            saveFig(f3, fname, FolderName);
            close(f3);
            
            %% Taken from Rhett's Code
            mvmt_good_starts = mvmt_pulses.starts(mvmt_good);
            rest_good_starts = rest_pulses.starts(rest_good);

            mvmt_good_mid = mvmt_good_starts + floor(mvmt_periods(mvmt_good)./2);
            rest_good_mid = rest_good_starts + floor(rest_periods(rest_good)./2);

            mvmt_good_cnst_range = [mvmt_good_mid-min_length./2;mvmt_good_mid+min_length./2];
            rest_good_cnst_range = [rest_good_mid-min_length./2;rest_good_mid+min_length./2];
            % rest periods
            volt_rest = [];
            hemo_rest = [];

            rest_epoch_num = length(rest_good_starts);
            move_epoch_num = length(mvmt_good_starts);

           if rest_epoch_num ~= 0
                for r = 1:rest_epoch_num
                    start = rest_good_cnst_range(1,r);
                    fin = rest_good_cnst_range(2,r);
                    
                    temp_volt = data_volt(:,:,start:fin);
                    temp_hemo = data_hemo(:,:,start:fin);
                    pupil_rest = pupil_data2x(start:fin);
                    
                    %% Get the Correlation maps for Voltage
                    % Rest period 
                    % Reshape data temp to 2d  
                    sXY = vsfp_data.out.sX * vsfp_data.out.sY;
                    data_rest_2d = reshape(temp_volt,[sXY,length(temp_volt)]);
                    % Correlation With Pupil At Baseline
                    corr_map_temp = fastCorr2(data_rest_2d', pupil_rest');
                    corr_map_reg_volt_rest = reshape(corr_map_temp,[100,100]);
                    maskImg = imresize(maskImg, [vsfp_data.out.sX, vsfp_data.out.sY]); 
                    maskImg = maskImg(:,:,1);
                    %apply mask on transformed data 
                    masked_corr_map_reg_volt_rest = imrotate(corr_map_reg_volt_rest, 270);
                    maskk = maskImg < 100;
                    masked_corr_map_reg_volt_rest(maskk)= 0;
                    
                    
                    %% Get the Correlation maps for Hemo
                    % Rest period 
                    % Reshape data temp to 2d  
                    data_rest_2d = reshape(temp_hemo,[sXY,length(temp_hemo)]);
                    % Correlation With Pupil At Baseline
                    corr_map_temp = fastCorr2(data_rest_2d', pupil_rest');
                    corr_map_reg_hemo_rest = reshape(corr_map_temp,[100,100]);
                    %apply mask on transformed data 
                    masked_corr_map_reg_hemo_rest = imrotate(corr_map_reg_hemo_rest, 270);
                    maskk = maskImg < 100;
                    masked_corr_map_reg_hemo_rest(maskk)= 0;
                
                    volt_rest = cat(3, volt_rest, masked_corr_map_reg_volt_rest);
                    hemo_rest = cat(3, hemo_rest, masked_corr_map_reg_hemo_rest);
                    
                    
                end 
                if move_epoch_num == 0
                    volt_mvmt = zeros(100,100);
                    hemo_mvmt = zeros(100,100);
                else
                    %movement periods
                    volt_mvmt = [];
                    hemo_mvmt = [];
                end
                for r = 1:move_epoch_num
                    start = mvmt_good_cnst_range(1,r);
                    fin = mvmt_good_cnst_range(2,r);
                    
                    temp_volt = squeeze(data_volt(:,:,start:fin));
                    temp_hemo = squeeze(data_hemo(:,:,start:fin));
                    pupil_mvmt = pupil_data2x(start:fin);

                    %% Corr for Voltage 
                    % Movement Period 
                    % Reshape data temp to 2d  
                    data_mvmt_2d = reshape(temp_volt,[sXY,length(temp_volt)]);
                    % Correlation With Pupil At Baseline
                    corr_map_temp = fastCorr2(data_mvmt_2d', pupil_mvmt');
                    corr_map_reg_volt_mvmt = reshape(corr_map_temp,[100,100]);
                    %apply mask on transformed data 
                    masked_corr_map_reg_volt_mvmt = imrotate(corr_map_reg_volt_mvmt, 270);
                    maskk = maskImg < 100;
                    masked_corr_map_reg_volt_mvmt(maskk)= 0;

                    % Movement Period Hemo
                    % Reshape data temp to 2d  
                    data_mvmt_2d = reshape(temp_hemo,[sXY,length(temp_hemo)]);
                    % Correlation With Pupil At Baseline
                    corr_map_temp = fastCorr2(data_mvmt_2d', pupil_mvmt');
                    corr_map_reg_hemo_mvmt = reshape(corr_map_temp,[100,100]);
                    %apply mask on transformed data 
                    masked_corr_map_reg_hemo_mvmt = imrotate(corr_map_reg_hemo_mvmt, 270);
                    maskk = maskImg < 100;
                    masked_corr_map_reg_hemo_mvmt(maskk)= 0;
                    
                    volt_mvmt = cat(3, volt_mvmt, masked_corr_map_reg_volt_mvmt);
                    hemo_mvmt = cat(3, hemo_mvmt, masked_corr_map_reg_hemo_mvmt);
                end 

                %% plot mean correlation map for volt rest and mvmt signal
                % get atlas to overlay
                fixed_map = imresize(fixed_map, [vsfp_data.out.sX, vsfp_data.out.sY]);
                M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts

                f1 = figure(1); 
                subplot(3,1,1)
                spec = squeeze(squeeze(mean(mean(temp_volt))));
                spectrogram(spec',[],[],[],Fs,'yaxis')
                title('Filtered GS Voltage')
                axis square
                
                subplot(3,1,2)
                imagesc(mean(volt_rest,3)), colorbar
                hold on
                a = imagesc(M);
                hold off
                alpha(a,.2)
                title('Rest Max Corr Coef Volt and Pupil')
                avgTrial_Volt_rest(:,:,k) = mean(volt_rest,3);
                axis square

                subplot(3,1,3)
                imagesc(mean(volt_mvmt,3)), colorbar
                hold on
                a = imagesc(M);
                hold off
                alpha(a,.2)
                title('Movement Max Corr Coef Volt and Pupil')
                avgTrial_Volt_mvmt(:,:,k) = mean(volt_mvmt,3);
                axis square

                %% plot correlation map for hemo rest and mvmt signal
                % get atlas to overlay

                f2 = figure(2); 
                subplot(3,1,1)
                spec = squeeze(squeeze(mean(mean(temp_volt))));
                spectrogram(spec',[],[],[],200,'yaxis')
                title('Filtered GS Hemo')
                axis square
                
                subplot(3,1,2)
                imagesc(mean(hemo_rest,3)), colorbar
                hold on
                a = imagesc(M);
                hold off
                alpha(a,.2)
                title('Rest Max Corr Coef Hemo and Pupil')
                axis square 
                avgTrial_Hemo_rest(:,:,k) = mean(hemo_rest,3);
                
                subplot(3,1,3)
                imagesc(mean(hemo_mvmt,3)), colorbar
                hold on
                a = imagesc(M);
                hold off
                alpha(a,.2)
                title('Movement Max Corr Coef Hemo and Pupil')
                axis square
                avgTrial_Hemo_mvmt(:,:,k) = mean(hemo_mvmt,3);

                % save figures
                fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_CorrMap_Hemo'];
                saveFig(f2, fname, FolderName);
                fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_CorrMap_Volt'];
                saveFig(f1, fname, FolderName);
                close all;
            end
        end
        
        mean_avgTrial_Hemo_mvmt = mean(avgTrial_Hemo_mvmt,3);
        mean_avgTrial_Hemo_rest = mean(avgTrial_Hemo_rest,3);
        
        mean_avgTrial_Volt_mvmt = mean(avgTrial_Volt_mvmt,3);
        mean_avgTrial_Volt_rest = mean(avgTrial_Volt_rest,3);
        
        % plot correlation map for avg volt signal
        f1 = figure(1); 
        subplot(2,1,1)
        imagesc(mean_avgTrial_Volt_rest), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Rest Max Corr Coef Volt and Pupil')
        axis square 
        
        subplot(2,1,2)
        imagesc(mean_avgTrial_Volt_mvmt), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Movement Max Corr Coef Volt and Pupil')
        axis square
        
        % plot correlation map for hemo signal
        f2 = figure(2); 
        subplot(2,1,1)
        imagesc(mean_avgTrial_Hemo_rest), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Rest Max Corr Coef Hemo and Pupil')
        axis square;
        
        subplot(2,1,2)
        imagesc(mean_avgTrial_Hemo_mvmt), colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title('Movement Max Corr Coef Hemo and Pupil')
        axis square

        fname = [mouse,'_','avg','_CorrMap_Hemo'];
        saveFig(f1, fname, FolderName);
        fname = [mouse,'_','avg','_CorrMap_Volt'];
        saveFig(f2, fname, FolderName);
        close all;
        
        %save all the avg corr maps for both hemo and voltage for each
        %animal    
        fname = [mouse,'_','avg','_CorrMap_Volt_rest'];
        save(fullfile(FolderName, fname), 'avgTrial_Volt_rest')
        fname = [mouse,'_','avg','_CorrMap_Volt_mvmt'];
        save(fullfile(FolderName, fname), 'avgTrial_Volt_mvmt')
        fname = [mouse,'_','avg','_CorrMap_Hemo_rest'];
        save(fullfile(FolderName, fname), 'avgTrial_Hemo_rest')
        fname = [mouse,'_','avg','_CorrMap_Hemo_mvmt'];
        save(fullfile(FolderName, fname), 'avgTrial_Hemo_mvmt')
        
    end   
end