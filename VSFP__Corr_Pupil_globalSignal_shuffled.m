% Function takes xlsx sheet and generates figures for following analysis
%       Modified from VSFP_Corr_Pupil_globalSignal.m 
%       
%       Takes 100 random pairs of imaging and pupil data to creat shuffled
%       control data
%      
% Inputs: 
%        optional mice number, or list of numbers i.e.
%        VSFP_Corr_Pupil_globalSignal('VSFP24','VSFP25') or VSFP_Corr_Pupil_globalSignal('VSFP24')
%        if no input is given, function runs through all animals
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP__Corr_Pupil_globalSignal_shuffled(varargin)
    close all
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

    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_Corr_Pupil_globalSignal_shuffled\'];   % Your destination folder


    % Get list of all good trials to pull pupil data from 
    all_mice = {'VSFP24', 'VSFP27', 'VSFP28', 'VSFP29', 'VSFP30'};
    all_trials_shuffle = [];
    all_dates_shuffle = [];
    for s = 1:length(all_mice)
        temp = find(contains(T.Mouse,all_mice{s}));
        all_trials_shuffle = [all_trials_shuffle;T.Trials(temp)];
        all_dates_shuffle = [all_dates_shuffle; T.Date(temp)];
    end

    lag_all = [];
    corr_all = [];
    lag_all_volt = [];
    corr_all_volt = [];
    lag_all_hemo = [];
    corr_all_hemo = [];

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

        % randomly pick another trial and load only the voltage imaging data from
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
        
        %Volt GS
        gsVolt = vsfp_data.out.gsVolt;
        clear vsfp_data

        % randomly pick another trial and load only the hemodynamic imaging data from
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
        
        % Hemo GS
        gsHemo = -1 * vsfp_data.out.gsHemo; %LP;
        
        % account for differences in length between trials by indexing based
        % on the shortest trace
        minL = min([length(pupil_data2x), length(gsVolt), length(gsHemo)]);
        pupil_data2x = pupil_data2x(1:minL);
        gsVolt = gsVolt(1:minL);
        gsHemo = gsHemo(1:minL);

        %% find the cross corr with global signal for the voltage data
        [r,lags] = xcorr(gsVolt,pupil_data2x,'normalized');
        Fs = vsfp_data.out.Fs;         % Sampling Frequency (Hz)
        t = lags/Fs;                   % Time Vector (seconds)
        
        
        % plot z-score of global voltage and pupil signal 
        f2 = figure(2);
        plot(gsVolt, 'DisplayName','Global Voltage Signal', 'LineWidth',2, 'Color', '#D95319')
        hold on
        plot(pupil_data2x, 'DisplayName','Pupil Diameter', 'LineWidth',2, 'Color', '#EDB120')
        xlabel( 'time (s)')
        ylabel( 'Z-Score')
        legend('Location','southeast')
        hold off 
        sgtitle('Correlation with Pupil Diameter (Volt/Pupil)')
        
        % plot corr between two signals in smaller window
        f3 = figure(3); 
        len = round(length(t)/ 2);
        ind = (len - 250: len + 250);
        x = t(ind);
        y = r(ind);
        plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
        hold on
        lag_all_volt = [lag_all_volt; x];
        corr_all_volt = [corr_all_volt;y];
        
        
        %% repeat but for hemo signal 
        [r,lags] = xcorr(gsHemo,pupil_data2x,'normalized');
        t = lags/Fs;                   % Time Vector (seconds)
        
        % plot z-score of global and pupil signal 
        f5 = figure(5);
        plot(gsHemo, 'DisplayName','Global Hemo Signal', 'LineWidth',2, 'Color', '#D95319')
        hold on
        plot(pupil_data2x, 'DisplayName','Pupil Diameter', 'LineWidth',2, 'Color', '#EDB120')
        xlabel( 'time (s)')
        ylabel( 'Z-Score')
        legend('Location','southeast')
        hold off 
        sgtitle('Correlation with Pupil Diameter (Hemo/Pupil)')
        
        % plot corr between two signals in smaller window
        f6 = figure(6); 
        len = round(length(t)/ 2);
        ind = (len - 250: len + 250);
        x = t(ind);
        y = r(ind);
        plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
        hold on
        lag_all_hemo = [lag_all_hemo; x];
        corr_all_hemo = [corr_all_hemo;y];

        %% correlation between volt and hemo
        [r,lags] = xcorr(gsHemo,gsVolt,'normalized');
        t = lags/Fs;   
        len = round(length(t)/ 2);
        ind = (len - 250: len + 250);
        x = t(ind);
        y = r(ind);
        f9 = figure(9);
        plot(x ,y, '--','color', [0.25, 0.25, 0.25], 'LineWidth',1)
        title('Correlation Between Hemo and Volt (Hemo/Volt)')
        hold on
        lag_all = [lag_all; x];
        corr_all = [corr_all;y];

        % save all the figures
        fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_zScore_globVolt_Pupil'];
        saveFig(f2, fname, FolderName);
        fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_zScore_globHemo_Pupil'];
        saveFig(f5, fname,FolderName);
        close (f2,f5);

        close all

    end
    %% plot the 100 trial averages
    % plot the cross corr for hemo/volt
    % add line to plot for avg signal across trials
    f9 = figure(9);
    plot(mean(lag_all,1), mean(corr_all,1), 'k', 'LineWidth',2)
    title( 'Trial Average Hemo and Volt');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')
    
    % calc standard error and plot that for each animal
    %err = std(corr_all,[],1)/sqrt(size(corr_all,1));
    err = std(corr_all,[],1)/sqrt(size(corr_all,1));
    f10 = figure(10);
    shade(mean(lag_all,1), mean(corr_all,1)+err,mean(lag_all,1),mean(corr_all,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all,1), mean(corr_all,1), 'k', 'LineWidth',2)
    ylim([-0.15 0.25])
    title( 'Trial Average Hemo and Volt');
    hold off
    
    %% plot for volt 
    % add line to plot for avg signal across trials
    f3 = figure(3);
    plot(mean(lag_all_volt,1), mean(corr_all_volt,1), 'k', 'LineWidth',2)
    title( 'Trial Average Volt and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')
    
    % calc standard error and plot that for each animal
    %err = std(corr_all,[],1)/sqrt(size(corr_all,1));
    err = std(corr_all_volt,[],1)/sqrt(size(corr_all_volt,1));
    f7 = figure(7);
    shade(mean(lag_all_volt,1), mean(corr_all_volt,1)+err,mean(lag_all_volt,1),mean(corr_all_volt,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all_volt,1), mean(corr_all_volt,1), 'k', 'LineWidth',2)
    ylim([-0.1 0.6])
    title( 'Trial Average Volt and Pupil');
    hold off
    
    %% plot for hemo 
    % add line to plot for avg signal across trials
    f6 = figure(6);
    plot(mean(lag_all_hemo,1), mean(corr_all_hemo,1), 'k', 'LineWidth',2)
    title( 'Trial Average Hemo and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')
    
    % calc standard error and plot that for each animal
    %err = std(corr_all,[],1)/sqrt(size(corr_all,1));
    err = std(corr_all_hemo,[],1)/sqrt(size(corr_all_hemo,1));
    f8 = figure(8);
    shade(mean(lag_all_hemo,1), mean(corr_all_hemo,1)+err,mean(lag_all_hemo,1),mean(corr_all_hemo,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all_hemo,1), mean(corr_all_hemo,1), 'k', 'LineWidth',2)
    ylim([-0.05 0.3])
    title( 'Trial Average Hemo and Pupil');
    hold off

    fname = 'Avg_standardError_Score_globVolt_Hemo';
    saveFig(f10, fname, FolderName); 
    fname = 'Avg_zScore_globVolt_Hemo';
    saveFig(f9,fname, FolderName); 
    
    fname = 'Avg_standardError_Score_globVolt_Pupil';
    saveFig(f7, fname, FolderName); 
    fname = 'Avg_zScore_globVolt_Pupil';
    saveFig(f3, fname, FolderName); 
    
    fname = 'Avg_standardError_Score_globHemo_Pupil';
    saveFig(f8, fname, FolderName); 
    fname = 'Avg_zScore_globHemo_Pupil';
    saveFig(f6, fname, FolderName); 
    close all;
end