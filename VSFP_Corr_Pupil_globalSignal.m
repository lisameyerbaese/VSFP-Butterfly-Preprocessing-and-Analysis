% Function takes xlsx sheet and generates figures for following analysis
%   1) calculates global signal and takes pupil diameter 
%      obtains correlation and time lag for both signals
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

function VSFP_Corr_Pupil_globalSignal(varargin)
    
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

    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_Corr_Pupil_globalSignal\'];   % Your destination folder

    lag_all_mice = [];
    corr_all_mice = [];
    lag_all_mice_volt = [];
    corr_all_mice_volt = [];
    lag_all_mice_hemo = [];
    corr_all_mice_hemo = [];
    fftVolt_all = [];
    fftHemo_all = [];
    
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        lag_all = [];
        corr_all = [];
        lag_all_volt = [];
        corr_all_volt = [];
        lag_all_hemo = [];
        corr_all_hemo = [];
        
        for k = 1:length(all_trials)

            %get corresponding pupil data
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;

            % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            % get pupil data
            pupil_data2x = zscore(vsfp_data.pupil_data2x);

            Fs = vsfp_data.out.Fs;         % Sampling Frequency (Hz)
            
            gsVolt = vsfp_data.out.gsVolt;
            gsHemo = -1 * vsfp_data.out.gsHemo; %LP;
            % low pass filter hemo signal at 1 Hz
            gsHemo = lowpass(gsHemo,1,Fs,ImpulseResponse="iir",Steepness=0.95);

            %gsVolt = calcGlobalSingal(vsfp_data.out.imgDR3 .* vsfp_data.out.mask);
            %gsHemo = calcGlobalSingal(-1 * vsfp_data.out.imgDS .* vsfp_data.out.mask);

            % load in the voltage and hemo data, and make it 2D, need this
            % for making the 2D plot of GS
            voltData2D = reshape(vsfp_data.out.imgDR3, 10000,[]);
            hemoData2D = reshape(-1 * vsfp_data.out.hemoLP, 10000, []);
            
            % find the cross corr with global signal for the voltage data
            [r,lags] = xcorr(gsVolt(1:length(pupil_data2x)),pupil_data2x,'normalized');
            t = lags/Fs;                   % Time Vector (seconds)
           
            % plot 2D image of global signal
            f1 = figure(1);
            imagesc(voltData2D), colorbar
            caxis([-0.05, 0.05])
            ylabel('Pixels from Ant to Post') 
            xlabel('time')
            sgtitle('Voltage Activity in Pixel Over Time')
            
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
            [r,lags] = xcorr(gsHemo(1:length(pupil_data2x)),pupil_data2x,'normalized');
            t = lags/Fs;                   % Time Vector (seconds)
            
            % plot 2D image of global signal
            f4 = figure(4);
            imagesc(hemoData2D), colorbar
            caxis([-0.05, 0.05])
            ylabel('Pixels from Ant to Post') 
            xlabel('time')
            title('Hemo Activity in Pixel Over Time')
            
            % plot z-score of global and pupil signal 
            f5 = figure(5);
            plot(gsHemo, 'DisplayName','Global Hemo Signal', 'LineWidth',2, 'Color', '#D95319')
            hold on
            plot(pupil_data2x, 'DisplayName','Pupil Diameter', 'LineWidth',2, 'Color', '#EDB120')
            xlabel( 'time (s)')
            ylabel( 'Z-Score')
            legend('Location','southeast')
            hold off 
            title('Correlation with Pupil Diameter (Hemo/Pupil)')
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_zScore_gsHemo_Pupil'];
            saveFig(f5, fname, FolderName);
            close(f5)

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
            
            %% plot the FFT for both the GS for Voltage and Hemo
            
            f9 = figure(9);
            Fs = vsfp_data.out.Fs;   % Sampling Frequency (Hz)    
            L = 4096;                % Use a fixed length to be able to avg across all GS
            f = Fs/L*(0:(L/2));
            
            % plot only from 0.7  Hz to 14Hz
            startInd = find(round(f,2) == 0.7);
            endInd = find(round(f,2) == 14);
            fftVolt = fft(gsVolt(1:L));
            P2 = abs(fftVolt/L);
            fftVolt = P2(1:L/2+1);
            fftVolt(2:end-1) = 2*fftVolt(2:end-1);
            subplot(2,1,1)
            plot(f(startInd(1):endInd(1)),fftVolt(startInd(1):endInd(1)),'k') 
            title('Single-Sided Amplitude Spectrum of GS Volt')
            xlabel("f (Hz)")
            ylabel("|P1(f)|")
            
            tempGSHemo = vsfp_data.out.gsHemo;
            fftHemo = fft(tempGSHemo(1:L));
            P2 = abs(fftHemo/L);
            fftHemo = P2(1:L/2+1);
            fftHemo(2:end-1) = 2*fftHemo(2:end-1);
            subplot(2,1,2)
            plot(f(startInd(1):endInd(1)),fftHemo(startInd(1):endInd(1)),'b') 
            title('Single-Sided Amplitude Spectrum of GS Hemo')
            xlabel("f (Hz)")
            ylabel("|P1(f)|")

            fftVolt_all = [fftVolt_all;fftVolt(startInd(1):endInd(1))];
            fftHemo_all = [fftHemo_all; fftHemo(startInd(1):endInd(1))];

            
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
            
        end
        lag_all_mice = [lag_all_mice; lag_all];
        corr_all_mice = [corr_all_mice; corr_all];
        
        lag_all_mice_volt = [lag_all_mice_volt; lag_all_volt];
        corr_all_mice_volt = [corr_all_mice_volt; corr_all_volt];
        
        lag_all_mice_hemo = [lag_all_mice_hemo; lag_all_hemo];
        corr_all_mice_hemo = [corr_all_mice_hemo; corr_all_hemo];
        
        %% plot for volt and hemo
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
        title( 'Trial Average Volt and Pupil');
        hold off
        
        %% plot for hemo 
        % add line to plot for avg signal across trials
        f6 = figure(f6);
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
        title( 'Trial Average Hemo and Pupil');
        hold off

        fname = [mouse,'_0', num2str(date) ,'_','avg','_standardError_Score_globVolt_Hemo'];
        saveFig(f10, fname, FolderName); 
        fname = [mouse,'_0', num2str(date) ,'_','avg','_zScore_globVolt_Hemo'];
        saveFig(f9,fname, FolderName); 
        
        %FolderName = 'X:\keilholz-lab\Lisa\VSFP ButterFly\Data\VSFP_Corr_Pupil_globalSignal\';   % Your destination folder
        %FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        fname = [mouse,'_0', num2str(date) ,'_','avg','_standardError_Score_globVolt_Pupil'];
        saveFig(f7, fname, FolderName); 
        fname = [mouse,'_0', num2str(date) ,'_','avg','_zScore_globVolt_Pupil'];
        saveFig(f3, fname, FolderName); 
        
        fname = [mouse,'_0', num2str(date) ,'_','avg','_standardError_Score_globHemo_Pupil'];
        saveFig(f8, fname, FolderName); 
        fname = [mouse,'_0', num2str(date) ,'_','avg','_zScore_globHemo_Pupil'];
        saveFig(f6, fname, FolderName); 
        close all;
    end 
    % avg for all animals for volt and hemo 
    f1 = figure(1);
    plot(mean(lag_all_mice,1), mean(corr_all_mice,1), 'k', 'LineWidth',2)
    title( 'Trial Average Hemo and Volt');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')

    f2 = figure(2);
    err = std(corr_all,[],1)/sqrt(size(corr_all,1));
    shade(mean(lag_all_mice,1), mean(corr_all_mice,1)+err,mean(lag_all_mice,1),mean(corr_all_mice,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all_mice,1), mean(corr_all_mice,1), 'k', 'LineWidth',2)
    title('Trial Average Hemo and Volt');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')
    
    fname = 'AvgAll_standardError_zScore_globVolt_Hemo';
    saveFig(f2, fname, FolderName);
    fname = 'AvgAll__zScore_globVolt_Hemo';
    saveFig(f1, fname, FolderName);
    close(f1, f2);
    
    
    % avg for all animals for volt 
    f1 = figure(1);
    plot(mean(lag_all_mice_volt,1), mean(corr_all_mice_volt,1), 'k', 'LineWidth',2)
    title( 'Trial Average Volt and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')

    f2 = figure(2);
    err = std(corr_all_volt,[],1)/sqrt(size(corr_all_volt,1));
    shade(mean(lag_all_mice_volt,1), mean(corr_all_mice_volt,1)+err,mean(lag_all_mice_volt,1),mean(corr_all_mice_volt,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all_mice_volt,1), mean(corr_all_mice_volt,1), 'k', 'LineWidth',2)
    title('Trial Average Volt and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')
    
    fname = 'AvgAll_standardError_zScore_globVolt_Pupil';
    saveFig(f2, fname, FolderName);
    fname = 'AvgAll_zScore_globVolt_Pupil';
    saveFig(f1, fname, FolderName);
    close (f1,f2);
    
    % avg for all animals for hemo 
    f1 = figure(1);
    plot(mean(lag_all_mice_hemo,1), mean(corr_all_mice_hemo,1), 'k', 'LineWidth',2)
    title( 'Trial Average Hemo and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')

    f2 = figure(2);
    err = std(corr_all_hemo,[],1)/sqrt(size(corr_all_hemo,1));
    shade(mean(lag_all_mice_hemo,1), mean(corr_all_mice_hemo,1)+err,mean(lag_all_mice_hemo,1),mean(corr_all_mice_hemo,1)- err, 'FillType',[1 2])
    hold on 
    plot(mean(lag_all_mice_hemo,1), mean(corr_all_mice_hemo,1), 'k', 'LineWidth',2)
    title( 'Trial Average Hemo and Pupil');
    ylabel('Correlation Coefficinet');
    xlabel('time (s)')

    % plot the avg FFT for the volt and hemo signal 
    f3 = figure(3);
    subplot(2,1,1)
    tempVolt = mean(fftVolt_all,1);
    tempHemo = mean(fftHemo_all,1);
    subplot(2,1,1)
    plot(f(59:1149),tempVolt,'k') 
    title('Single-Sided Amplitude Spectrum of GS Volt')
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    subplot(2,1,2)
    plot(f(59:1149),tempHemo,'b') 
    title('Single-Sided Amplitude Spectrum of GS Hemo')
    xlabel("f (Hz)")
    ylabel("|P1(f)|")
    
    fname = 'AvgAll_FFT';
    saveFig(f3, fname, FolderName);
    fname = 'AvgAll_standardError_zScore_globHemo_Pupil';
    saveFig(f2, fname, FolderName);
    fname = 'AvgAll_zScore_globHemo_Pupil';
    saveFig(f1, fname, FolderName);

    close all;
end

function [globalSignal] = calcGlobalSingal(data3D)
    [xDim, yDim, zDim] = size(data3D);
    data3D(isnan(data3D))=0;
    data2D = reshape(data3D, xDim * yDim, zDim);
    data2D = zscore(data2D, 0, 2);
    globalSignal = nanmean(data2D, 1);
end