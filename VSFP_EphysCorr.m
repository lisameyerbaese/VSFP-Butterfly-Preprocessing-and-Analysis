% Function loads in the preproccessed VSFP ephys data and plots the
% correlation between the LFP data and both the voltage and hemodynamic
% responses
%
% Written by Lisa Meyer-Baese

function VSFP_EphysCorr()

    close all;
    %get root of file name based on computer being used 
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    folderName = [startFile,'\VSFP ButterFly\Data\VSFP_EphysCorr'];
    %check to see if Folder Exists to Save data, if not create it
     if ~exist(folderName, 'dir')
       mkdir(folderName)
    end

    allTrials = dir([startFile,'\VSFP ButterFly\Data\VSFP_Ephys\']);


    % initalize variables for trail avg 
    all_CorrVolt = [];
    all_corrGSVolt = [];
    all_corrROIVolt = [];
    all_lagGSVolt = [];
    all_lagROIVolt = [];
    all_CorrHemo = [];
    all_corrGSHemo = [];
    all_corrROIHemo = [];
    all_lagGSHemo = [];
    all_lagROIHemo = [];
    all_CrossCorrGSHemoVolt = [];
    all_CrossCorrGSLFPVolt = [];
    all_CrossCorrGSLFPHemo = [];
    all_CrossCorrROILFPVolt = [];
    all_CrossCorrROILFPHemo = [];

    % load in trials of comparable depth between VSFPE3 and VSFPE4 
    %   VSFPE3, Trials 311 - 329, depth 1152 um
    %   VSFPE4, Trials 261 - 270, depth 1150 um
    %   trials = 3:20,86:95;
    %   trials = 8:20,86:95;

    % or load in trials only from one animal for a single depth 
    trials = 21:30; % VSFPE3, trials 361-370 with depth 1500um 
    for t = 1:length(trials)
        i = trials(t);
        %i = 29;
        load([allTrials(i).folder,'\',allTrials(i).name])
        fileName = allTrials(i).name(1:end-4);

        fs = out.Fs;
        L = length(out.gsVolt);
        T = 1/fs;
        time = (0:L-1)*T;     % Time vector

        % get a correlation map between imaging data and LFP data
        voltData = out.imgDR3; %imrotate(out.imgDR3, 270);
        hemoData = -1 * out.imgDS;
        [avgCorrVolt, corrGSVolt, corrROIVolt, lagGSVolt, lagROIVolt, crossCorrROIVolt, crossCorrGSVolt, avg_ROIVolt] = makeMap(voltData, out, sortedLFP, time, [fileName,'_Volt'], folderName, 'Volt'); 
        [avgCorrHemo, corrGSHemo, corrROIHemo, lagGSHemo, lagROIHemo, crossCorrROIHemo, crossCorrGSHemo, avg_ROIHemo] = makeMap(hemoData, out, sortedLFP, time, [fileName,'_Hemo'], folderName, 'Hemo'); 

        % plot cross correlation between global volt and hemo, GS volt/LFP
        % and hemo/lfp and ROI volt/lfp and ROI hemo/lfp
        
        % low pass filter voltGS at 0.1Hz
        voltGS = out.gsVolt;

        hemoGS = -1 * out.gsHemoLP; % invert the hemodynamic signal so that + = more HbT

        % low pass filter hemoGS at 1Hz
        hemoGS = lowpass(hemoGS,0.1,fs,ImpulseResponse="iir",Steepness=0.95);
        [cGSHemoVolt,lags] = xcorr(hemoGS,voltGS, (5*fs), 'normalize');

        f1 = figure(1);
        subplot(2,1,1)
        plot(time, zscore(voltGS), 'k','DisplayName','Volt')
        hold on 
        plot(time, zscore(hemoGS),'g', 'DisplayName','Hemo')
        xlabel('time (s)')
        legend()
        subplot(2,1,2)
        plot(lags/fs, cGSHemoVolt)
        title('Cross Corr (hemo/volt)')
        sgtitle(strrep( fileName , '_' , ' '))
        saveFig(f1, ['CrossCorr_VoltHemo_',fileName], folderName)
        close(f1)

        f2 = figure(2);
        subplot(2,2,1)
        plot(lags/fs, crossCorrROIVolt)
        title('Cross Corr ROI LFP/Volt')
        subplot(2,2,2)
        plot(lags/fs, crossCorrGSVolt)
        title('Cross Corr GS LFP/Volt')
        subplot(2,2,3)
        plot(lags/fs, crossCorrROIHemo)
        title('Cross Corr ROI LFP/Hemo')
        subplot(2,2,4)
        plot(lags/fs, crossCorrGSHemo)
        title('Cross Corr GS LFP/Hemo')
        sgtitle(strrep( fileName , '_' , ' '))
        saveFig(f2, ['CrossCorr_with_LFP_',fileName], folderName)
        close(f2) 

        f3 = figure(3);
        subplot(4,1,1)
        plot(time,avg_ROIVolt)
        subplot(4,1,2)
        plot(time(1:200*5), avg_ROIVolt(1:200*5))
        subplot(4,1,3) 
        plot(time, avg_ROIHemo)
        subplot(4,1,4) 
        plot(time(1:200*5), avg_ROIHemo(1:200*5))

        % save the values for each trial
        % first the cross corr values
        all_CrossCorrGSHemoVolt = [all_CrossCorrGSHemoVolt; cGSHemoVolt];
        all_CrossCorrGSLFPVolt = [all_CrossCorrGSLFPVolt;  crossCorrGSVolt];
        all_CrossCorrGSLFPHemo = [all_CrossCorrGSLFPHemo; crossCorrGSHemo];
        all_CrossCorrROILFPVolt = [all_CrossCorrROILFPVolt;  crossCorrROIVolt];
        all_CrossCorrROILFPHemo = [all_CrossCorrROILFPHemo; crossCorrROIHemo];

        % then the voltage data
        all_CorrVolt = cat(4, all_CorrVolt, avgCorrVolt);
        all_corrGSVolt = [all_corrGSVolt; corrGSVolt];
        all_corrROIVolt = [all_corrROIVolt; corrROIVolt];
        all_lagGSVolt = [ all_lagGSVolt;lagGSVolt ];
        all_lagROIVolt = [all_lagROIVolt; lagROIVolt ];
        
        % then the hemo
        all_CorrHemo = cat(4, all_CorrHemo, avgCorrHemo);
        all_corrGSHemo = [all_corrGSHemo;corrGSHemo ];
        all_corrROIHemo = [all_corrROIHemo;corrROIHemo ];
        all_lagGSHemo = [ all_lagGSHemo; lagGSHemo ];
        all_lagROIHemo = [all_lagROIHemo;lagROIHemo ];
   
    end

    % look at the average activity across all trials
    plotAvg(all_CorrVolt, all_corrGSVolt, all_corrROIVolt, all_lagGSVolt, all_lagROIVolt, '_Volt', folderName)
    plotAvg(all_CorrHemo, all_corrGSHemo, all_corrROIHemo, all_lagGSHemo, all_lagROIHemo, '_Hemo', folderName)

    % plot the average cross correlation between with the standard error 
    f2 = figure(2);
    subplotTotal = 5;
    plotAvgCrossCorr(lags, fs, all_CrossCorrGSHemoVolt, subplotTotal, 1, 'GS Hemo', 'Volt')
    plotAvgCrossCorr(lags, fs, all_CrossCorrGSLFPVolt, subplotTotal, 2, 'LFP', 'GS Volt')
    plotAvgCrossCorr(lags, fs, all_CrossCorrGSLFPHemo, subplotTotal, 3, 'LFP', 'GS Hemo')
    plotAvgCrossCorr(lags, fs, all_CrossCorrROILFPVolt, subplotTotal, 4, 'LFP', 'ROI Volt')
    plotAvgCrossCorr(lags, fs, all_CrossCorrROILFPHemo, subplotTotal, 5, 'LFP', 'ROI Hemo')
    saveFig(f2, 'AvgCrossCorr_', folderName)

    f3 = figure(3);
    plotAvgCrossCorr(lags, fs, all_CrossCorrGSHemoVolt, 1, 1, 'GS Hemo', 'Volt')
    saveFig(f3, 'AvgGS_Volt_Hemo_CrossCorr_', folderName)

    save([folderName,'\AvgVoltData.mat'],'all_CorrVolt', 'all_corrGSVolt', 'all_corrROIVolt', 'all_lagGSVolt', 'all_lagROIVolt')
    save([folderName,'\AvgHemoData.mat'], 'all_CorrHemo', 'all_corrGSHemo', 'all_corrROIHemo', 'all_lagGSHemo', 'all_lagROIHemo')

end

function plotAvg(total_Corr, corrGS, corrROI, lagGS, lagROI, fileName, folderName)

    % this function takes avg activity and makes plots for the correlation with activity from 4 Channels   
    %[~,trials] = size(corrROI);
    % select image used for mask
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));

    % load in cortical map Cortical_Map.png
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));

    % plot 4 individual Correlation maps for the 4 Channels of LFP
    % data
    ch = [8,16,24,32];
    colors = {'#0072BD','#D95319', '#EDB120', '#7E2F8E'};
    for i = 1:length(ch)

        % plot avg spatial map
        f5 = figure(5);
        maskImg = imresize(maskImg, [100, 100]); 
        maskImg = maskImg(:,:,1);
        %apply mask on transformed data 
        maskk = maskImg < 100;
        % avg across trials 
        temp = mean(total_Corr,4);
        temp = temp(:,:,ch(i));
        temp = imrotate(temp, 270);
        temp = imgaussfilt(temp,0.5);
        temp(maskk)= 0;
        fixed_map = imresize(fixed_map, [100,100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
        imagesc(temp), axis off square
        colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title([' Avg Map for Ch',num2str(ch(i))])

        saveFig(f5,['AvgMap_Ch',num2str(ch(i)),'_',fileName], folderName )
        close(f5)

        % plot scatter plot with all point and avg for the Correlation
        % values 
        f6 = figure(6);
        avgCorrGS = mean(corrGS,1);
        avgLagsGS = mean(lagGS,1);
        subplot(2,1,1)
        scatter(ch, corrGS(:, ch),'.', 'MarkerFaceColor',colors{i})
        hold on 
        scatter(ch(i), avgCorrGS(ch(i)),'k', "o", 'filled')
        ylabel('Max Cross Corr Value')
        xlabel('LFP Channel')
    
        subplot(2,1,2)
        scatter(ch(i), lagGS(:, ch(i)),".", 'MarkerFaceColor',colors{i})
        hold on 
        scatter(ch(i), avgLagsGS(ch(i)),'k', "o", 'filled')
        ylabel('Lag Value (s)')
        xlabel('LFP Channel')
        sgtitle(['Global Signal ',fileName])

        f7 = figure(7);
        avgCorrROI = mean(corrROI,1);
        avgLagsROI = mean(lagROI,1);
        subplot(2,1,1)
        scatter(ch, corrROI(:, ch),'.', 'MarkerFaceColor',colors{i})
        hold on 
        scatter(ch(i), avgCorrROI(ch(i)),'k', "o", 'filled')
        ylabel('Max Cross Corr Value')
        xlabel('LFP Channel')
    
        subplot(2,1,2)
        scatter(ch(i), lagROI(:, ch(i)),".", 'MarkerFaceColor',colors{i})
        hold on 
        scatter(ch(i), avgLagsROI(ch(i)),'k', "o", 'filled')
        ylabel('Lag Value (s)')
        xlabel('LFP Channel')
        sgtitle(['ROI Activity ',fileName])
    end

    saveFig(f6,['Scatter_CorrLags_GS_',fileName], folderName)
    saveFig(f7,['Scatter_CorrLags_ROI_',fileName], folderName)
    close(f6, f7)
end

function [avgCorr, corr_GS, corr_ROI, lag_GS, lag_ROI, crossCorr_ROI, crossCorr_GS, avg_ROI]  = makeMap(imgData, out,sortedLFP, time, fileName, FolderName, type)

    % select image used for mask
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));

    % load in cortical map Cortical_Map.png
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));

    % process imaging Data
    data2D = reshape(imgData,[10000,4096]);
    if isequal(type, 'Volt')
        img_bpfilt = make_ChebII_filter(1, out.Fs, [0.5,50], [0.5*0.9 50*1.1], 20);
        vsfp_data_filt =  filtfilt(img_bpfilt.numer,img_bpfilt.denom,data2D')';
    else 
         vsfp_data_filt =  data2D;
    end
    %vsfp_data_filt = data2D;
    avg_global = mean(vsfp_data_filt,1);
    avg_global = zscore(avg_global,[],2);
    %process ROI data
    % take a 2 by 2 pixel ROI [93:94,50:51]
    mask_new = zeros(100,100);
%     mask_new(93,51) = 1;
%     mask_new(94,51) = 1;
%     mask_new(93,52) = 1;
%     mask_new(94,52) = 1;
    mask_new(7,51) = 1;
    mask_new(8,51) = 1;
    mask_new(7,52) = 1;
    mask_new(8,52) = 1;
    ROI_data = imgData .* mask_new;
    vsfp_ROI_2d = reshape(ROI_data,[10000,4096]);
    if isequal(type, 'Volt')
        vsfp_ROI_filt =  filtfilt(img_bpfilt.numer,img_bpfilt.denom,vsfp_ROI_2d')';
    else 
        vsfp_ROI_filt = vsfp_ROI_2d;
    end
    %vsfp_ROI_filt = vsfp_ROI_2d;
    avg_ROI = mean(vsfp_ROI_filt,1);
    avg_ROI = zscore(avg_ROI,[],2);
    corr_GS = [];
    corr_ROI = [];
    lag_GS = [];
    lag_ROI = [];
    crossCorr_ROI = [];
    crossCorr_GS = [];
    
%     % process LFP data
%     if isequal(type, 'Volt')
%         % Ephys Low Pass Filter at 50 Hz
%         ephys_lpfilt = make_ChebII_filter(2, 20000, 50, 50*1.1, 20);
%         %ephys_lpfilt = make_ChebII_filter(1, 20000, [0.5,40], [0.5*0.9 40*1.1], 20);
%         sortedLFP = filtfilt(ephys_lpfilt.numer, ephys_lpfilt.denom,sortedLFP);
%         LFP_gs_data = sortedLFP';
%         LFP_ds = downsample(LFP_gs_data',100);
%         LFP_ds = -1 * LFP_ds'; %flip over 
% 
%     % if signal is Hemo, add some extra pre-processing steps to the LFP
%     % signal before calculating the correlation 
%     elseif isequal(type, 'Hemo')
%         % Ephys Low Pass Filter at 5Hz
%         ephys_lpfilt = make_ChebII_filter(2, 20000, 5, 50*1.1, 20);
%         sortedLFP = filtfilt(ephys_lpfilt.numer, ephys_lpfilt.denom,sortedLFP);
%         LFP_gs_data = sortedLFP';
%         LFP_ds = downsample(LFP_gs_data',100);
%         LFP_ds = -1 * LFP_ds'; %flip over 
%     end
    
    % Ephys Low Pass Filter at 50 Hz
    ephys_lpfilt = make_ChebII_filter(2, 20000, 50, 50*1.1, 20);
    %ephys_lpfilt = make_ChebII_filter(1, 20000, [0.5,40], [0.5*0.9 40*1.1], 20);
    sortedLFP = filtfilt(ephys_lpfilt.numer, ephys_lpfilt.denom,sortedLFP);
    LFP_gs_data = sortedLFP';
    LFP_ds = downsample(LFP_gs_data',100);
    LFP_ds = -1 * LFP_ds'; %flip over 

    % Faster correlation functions
    for i2 = 1:32
        % Correlation between voltage and LFP 
        corr_map_temp = fastCorr2(vsfp_data_filt',(LFP_ds(i2,:))');
        corr_map_ephys(:,:,i2) = reshape(corr_map_temp,[100,100]);
        
        % Correlation between GS and LFP, look at max correlation within a 5s window 
        [cGS,lags] = xcorr(LFP_ds(i2,:),avg_global,'normalized', out.Fs * 5);
        [val, ind] = max(cGS);
        corr_GS(i2) = val;
        lag_GS(i2) = lags(ind)/out.Fs;

        % Correlation between ROI and LFP 
        [cROI,lags] = xcorr(LFP_ds(i2,:),avg_ROI,'normalized', out.Fs * 5);
        [val, ind] = max(cROI);
        corr_ROI(i2) = val; 
        lag_ROI(i2) = lags(ind)/out.Fs;

        if isequal(i2, 8) % if equal to channel 8, save the corr traces
            crossCorr_ROI = [crossCorr_ROI; cROI];
            crossCorr_GS = [crossCorr_GS; cGS];
        end     
    end
    
    % Display Corr Map for GS
    f2 = figure(2);
    imagesc(corr_GS), colorbar, caxis([-0.5 0.5])
    axis off
    sgtitle(strrep( fileName , '_' , ' '))
    sgtitle(['Corr with GS ', strrep( fileName , '_' , ' ')])
  
    % Display Corr Map for ROI
    f3 = figure(3);
    imagesc(corr_ROI), colorbar, caxis([-0.5 0.5])
    axis off
    sgtitle(['Corr with ROI ', strrep( fileName , '_' , ' ')])
    
    %Plot Corr_Map_Figure
    f4 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i2 = 1:32
        maskImg = imresize(maskImg, [100, 100]); 
        maskImg = maskImg(:,:,1);
        %apply mask on transformed data 
        maskk = maskImg < 100;
        temp = imrotate(corr_map_ephys(:,:,i2),270);
        temp = imgaussfilt(temp,0.5);
        temp(maskk)= 0;
        fixed_map = imresize(fixed_map, [100,100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
        
        subplot(4,8,i2)
        imagesc(temp), axis off square
        colorbar
        hold on
        a = imagesc(M);
        hold off
        alpha(a,.2)
        title(['Ch',num2str(i2)])
        
    end
    sgtitle(strrep( fileName , '_' , ' '))
    
    %Plot LFP Traces with GS and ROI Signal
    f5 = figure(5);
    plot(time, zscore(avg_global) + (36 *2), 'r', 'DisplayName','GS')
    hold on 
    plot(time, zscore(avg_ROI) + (33 *2), 'b','DisplayName','ROI Signal')
    legend()
    hold on 

    plot(time, zscore(LFP_ds(8,:)) + ((32-8) *2),'color', 'k', 'DisplayName','Ch8')
    plot(time, zscore(LFP_ds(16,:)) + ((32-16) *2),'color', 'k', 'DisplayName','Ch16')
    plot(time, zscore(LFP_ds(24,:)) + ((32-24) *2), 'color','k', 'DisplayName','Ch34')
    plot(time, zscore(LFP_ds(32,:)) + ((32-32) *2), 'color','k', 'DisplayName','Ch32')
    hold on
    sgtitle(strrep( fileName , '_' , ' '))
    avgCorr = corr_map_ephys;  %return the correlation map for each trial 

    %save all the figures made, and then close them 
    saveFig(f2, ['CorrMap_GS_',fileName], FolderName)
    saveFig(f3, ['CorrMap_ROI_',fileName], FolderName)
    saveFig(f4, ['Spatial_CorrMap_ROI_', fileName], FolderName)
    saveFig(f5, ['Raw_Traces_',fileName], FolderName)
    close(f2,f3,f4,f5)


end

function plotAvgCrossCorr(lags, fs, avgCrossCorr, subplotTotal, subplotInd, firstTrace, secondTrace)
    % this function takes the cross corr values and plots the avg along
    % with the standard error bars
    subplot(subplotTotal, 1, subplotInd)
    SE = std(avgCrossCorr,[],1)/sqrt(size(avgCrossCorr,1));
    shade(lags/fs, mean(avgCrossCorr,1) + SE,lags/fs,mean(avgCrossCorr,1) - SE, 'FillType',[1 2])
    hold on
    plot(lags/fs,mean(avgCrossCorr,1))
    [val, ind] = max(mean(avgCrossCorr,1));
    xline((lags(ind))/fs);
    xlabel('time (s)')
    title(['Cross Correlation between ', firstTrace,' and ', secondTrace,' Max:', num2str(val), ' Time:', num2str((lags(ind))/fs)])
end
