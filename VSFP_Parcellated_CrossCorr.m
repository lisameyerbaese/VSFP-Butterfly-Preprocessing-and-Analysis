% Function takes xlsx sheet and generates figures for following analysis
%   1) Takes seperate hemo and voltage signal, and averages the activity within parcells of the 
%      Allen Brain Atlas, then finds the cross corr for each region 
%   2) Edited from an older version called VSFP_Delay_CorrCorr.m
%  
% Inputs: 
%       method = type of data to use, either projection or ratiometric
%       type = whether or not you want to do GSR on the hemo data prior to
%              calculating the cross corr 
%
%       VSFP_Parcellated_CrossCorr('projection', 'yGSR')
%       VSFP_Parcellated_CrossCorr('ratiometric', 'nGSR')
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP_Parcellated_CrossCorr(method, type)
    
    close all;
    %get root of file name based on computer being used 
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile,'\VSFP ButterFly\Info'])
    addpath ([startFile,'\VSFP ButterFly\Code\Analysis\Lake Code'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    FolderName = [startFile,'\VSFP ButterFly\Data\VSFP_Parcellated_CrossCorr\',method,'_',type,'\'];   % Your destination folder for saving figures
    %check to see if Folder Exists to Save data, if not create it
     if ~exist(FolderName, 'dir')
       mkdir(FolderName)
    end
    
    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % delete VSFP 25 no pupil video


    % load in masked ROIs 
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas\HRF ROIs';
    ROIs = dir(path_map);

    % load in whole cortical map
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    
    disp('Selecting Cortical Map to Use for Registration')
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));
    
    % initalize empty parameters for avg across animals
    allTraces = [];
    allCorr = [];
    allCorrVP = [];
    allCorrHP = [];
    allLags = [];
    maxValVH = {};
    maxIndVH = {};
    maxValVP = {};
    maxIndVP = {};
    maxValHP = {};
    maxIndHP = {};
    maxLags = [];
    maxCorr = [];

    miceColor = {'r','b','g','y','m'};

    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        % initalize empty parameters for trial avg data
        avgCorr = [];
        avgCorrVP = [];
        avgCorrHP = [];
        avgLags = [];
        
        for k = 1:length(all_trials)
        
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;

             % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            % get pupil data
            pupilData = zscore(vsfp_data.pupil_data2x); 

            % load in the voltage and hemo data, and rotate it so that
            % midline is vertial not horizontal 
            if isequal(method, 'projection')
                %run on the projection method data
                volt = imrotate(vsfp_data.out.projectedVolt, 270);
                hemo = imrotate((-1 * vsfp_data.out.projectedHemo), 270);
                volt(isnan(volt))=0;
                hemo(isnan(hemo))=0;
                % low pass filter the hemo signal to 5Hz 
                hemo2D = reshape(hemo, 10000,[]); 
                lpHemo2D = lowpass(hemo2D',5,vsfp_data.out.Fs); %filters each column independantly
                hemo = reshape(lpHemo2D',100,100,[]);
                tempMask = imrotate(vsfp_data.out.mask, 270);
                gsHemo = calcGlobalSingal(hemo.* tempMask);
            else
                volt = imrotate(vsfp_data.out.imgDR3, 270);
                hemo = imrotate((-1 * vsfp_data.out.hemoLP), 270);
                gsHemo = -1 * vsfp_data.out.gsHemoLP;
            end
         
            % generate time vector 
            Fs = vsfp_data.out.Fs;     % Sampling Frequency (Hz)
            L = size(volt,3);
            t = (0:L-1)/Fs;        % Time vector

            roi = 1;
            crossCorr = [];
            crossCorrVP = [];
            crossCorrHP = [];
            lags = [];

            % loop through each ROI/parcell 
            for r = 3:length(ROIs)
                
                %load in the mask 
                tempMask = imread(fullfile(ROIs(r).folder, ROIs(r).name));
                tempMask = imresize(tempMask(:,:,3),[100,100]);

                maskk = tempMask <= 60;
                masked_volt = volt .* maskk;
                masked_hemo = hemo .* maskk;
                
                masked_volt2D = mean(reshape(masked_volt, 10000,[]),1);
                masked_hemo2D = mean(reshape(masked_hemo, 10000,[]),1);
                
                % perform global signal regression on hemo data if selected
                % as input
                if isequal(type, 'yGSR')
                    [masked_hemo2D] = GSR(masked_hemo2D, gsHemo);
                end

                x_f = masked_volt2D-mean(masked_volt2D);     % normalize x
                x_f = x_f/norm(x_f);
                y_f = masked_hemo2D-mean(masked_hemo2D);
                y_f = y_f/norm(y_f);                         % normalize y

                % add a 30 sec window to smooth both traces
%                 window = 3 * Fs;
%                 x_f = zscore(movmean(gsVoltTemp,window)); 
%                 y_f = zscore(movmean(gsHemoTemp,window));

                %plot the two signals 
                f4 = figure(4);
                subplot(3,3,roi)
                plot(t,x_f)
                hold on 
                plot(t,y_f)
                %plot(t(1:length(pupilData)),  pupilData, 'k')
                title(ROIs(r).name(1:end-4))
                xlabel('time (s)')
                ylabel('norm dF/F')
                roi = roi +1;
                sgtitle('Avg Trace per ROI: Volt(blue)  Hemo(orange)')

                %look also at the cross correlation between the two signals
                % going out +/- 5 seconds
                [corrVal,lagVal] = xcorr(y_f, x_f, 250, 'normalize');
                [corrValVP,lagValVP] = xcorr(x_f(1:length(pupilData)),pupilData, 250, 'normalize');
                [corrValHP,lagValHP] = xcorr(y_f(1:length(pupilData)),pupilData, 250, 'normalize');
                
                crossCorr = [crossCorr; corrVal];
                crossCorrVP = [crossCorrVP; corrValVP];
                crossCorrHP = [crossCorrHP; corrValHP];
                lags = [lags; lagVal/Fs];
               
            end

            % Cross Corr traces
            f5 = figure('units','normalized','outerposition',[0 0 1 1]);
            for i = 1:9
                subplot(2,5,i)
                [maxVal,maxInd] = max(abs(crossCorr(i,:)));
                plot(lags(i,:),crossCorr(i,:),'DisplayName', 'Voltage and Hemo')
                hold on
                temp = lags(i,:);
                xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

                [maxVal,maxInd] = max(abs(crossCorrVP(i,:)));
                plot(lags(i,:),crossCorrVP(i,:), 'DisplayName', 'Voltage and Pupil')
                hold on
                temp = lags(i,:);
                xline(temp(maxInd),'r','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

                [maxVal,maxInd] = max(abs(crossCorrHP(i,:)));
                plot(lags(i,:),crossCorrHP(i,:), 'DisplayName', 'Hemo and Pupil')
                hold on
                temp = lags(i,:);
                xline(temp(maxInd),'y','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

                legend('Location','northoutside')
                axis square
                xlabel('time (s)')
                ylabel('Correlation')
                title([ROIs(i+2).name(1:end-4)])

            end

            % Cross Corr traces in only the postitive lag window
            f6 = figure('units','normalized','outerposition',[0 0 1 1]);
            for i = 1:9
                subplot(2,5,i)
                tempInd = find(lags(i,:) == 0);
                startInd = tempInd - 250;
                endInd = tempInd + 250;
                [maxVal,maxInd] = max(crossCorr(i,startInd:endInd));
                plot(lags(i,:),crossCorr(i,:),'DisplayName', 'Voltage and Hemo')
                hold on
                temp = lags(i,startInd:endInd);
                xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

                legend('Location','northoutside')
                axis square
                xlabel('time (s)')
                ylabel('Correlation')
                title([ROIs(i+2).name(1:end-4)])

            end
            sgtitle('Cross Corr Between Voltage and Hemo showing max within +/- 5s')

            %save avg figure for each trial
            fname = [mouse,'_0',num2str(date),'_',num2str(trial),'_ROI_Traces'];
            saveFig(f4, fname, FolderName)
            fname = [mouse,'_0',num2str(date),'_',num2str(trial),'_Cross_Corr_Traces'];
            saveFig(f5, fname, FolderName)
            fname = [mouse,'_0',num2str(date),'_',num2str(trial),'_Cross_Corr_VoltHemo'];
            saveFig(f6, fname, FolderName)
            close (f4, f5, f6);

            % save individual trial to do a trial avg
            avgCorr = cat(3, avgCorr, crossCorr);
            avgCorrHP = cat(3, avgCorrHP, crossCorrHP);
            avgCorrVP = cat(3, avgCorrVP, crossCorrVP);
            avgLags = cat(3,avgLags, lags);
        end
        
        % Plot Animal Specific Avg Data
        % Cross Corr
        f4 = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:9
            subplot(2,5,i)
            tempCorr = mean(avgCorr(i,:,:),3);
            [maxVal,maxInd] = max(abs(tempCorr));
            plot(lags(i,:),tempCorr, 'DisplayName', 'Voltage and Hemo')
            hold on
            temp = lags(i,:);
            xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

            tempCorr = mean(avgCorrVP(i,:,:),3);
            [maxVal,maxInd] = max(abs(tempCorr));
            plot(lags(i,:),tempCorr, 'DisplayName', 'Voltage and Pupil')
            hold on
            temp = lags(i,:);
            xline(temp(maxInd),'r','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

            tempCorr = mean(avgCorrHP(i,:,:),3);
            [maxVal,maxInd] = max(abs(tempCorr));
            plot(lags(i,:),tempCorr, 'DisplayName', 'Hemo and Pupil')
            hold on
            temp = lags(i,:);
            xline(temp(maxInd),'y','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

            legend('Location','northoutside')
            axis square
            xlabel('time (s)')
            ylabel('Correlation')
            title([ROIs(i+2).name(1:end-4)])
        end

        f8 = figure(8);
        for i = 1:9
            [~,~,trials] = size(avgCorr);
            subplot(2,5,i)
            axis square
            hold on
                for j = 1:trials
                    tempCorr = avgCorr(i,:,j); 
                    [~,maxInd] = max(abs(tempCorr));
                    temp = lags(i,:);
                    scatter(temp(maxInd),tempCorr(maxInd),'filled',miceColor{m})
                    maxLags = [maxLags, temp(maxInd)];
                    maxCorr = [maxCorr, tempCorr(maxInd)];
                    hold on
                    xlabel('Lag (s)')
                    ylabel('Corr')
                end
                title([ROIs(i+2).name(1:end-4)])
        end

        f9 = figure(9);
        for i = 1:9
            [~,~,trials] = size(avgCorrVP);
            subplot(2,5,i)
            axis square
            hold on
                for j = 1:trials
                    tempCorr = avgCorrVP(i,:,j);
                    [~,maxInd] = max(abs(tempCorr));
                    temp = lags(i,:);
                    scatter(temp(maxInd),tempCorr(maxInd),'filled',miceColor{m})
                    hold on
                    xlabel('Lag (s)')
                    ylabel('Corr')
                end
            title([ROIs(i+2).name(1:end-4)])
        end

        f10 = figure(10);
        for i = 1:9
            [~,~,trials] = size(avgCorrHP);
            subplot(2,5,i)
            axis square
            hold on
                for j = 1:trials
                    tempCorr = avgCorrHP(i,:,j);
                    [~,maxInd] = max(abs(tempCorr));
                    temp = lags(i,:);
                    scatter(temp(maxInd),tempCorr(maxInd),'filled',miceColor{m})
                    hold on
                    xlabel('Lag (s)')
                    ylabel('Corr')
                end
            title([ROIs(i+2).name(1:end-4)])
        end

        f11 = figure('units','normalized','outerposition',[0 0 1 1]);
        for i = 1:9

            subplot(2,5,i)
            %tempInd = find(lags(i,:) == 0); % look only at the positive lag window
            tempInd = find(lags(i,:) == 0);
            startInd = tempInd - 250;
            endInd = tempInd + 250;
            tempCorr = mean(avgCorr(i,:,:),3);
            [maxVal,maxInd] = min(tempCorr(tempInd:end));
            plot(lags(i,tempInd:end),tempCorr(tempInd:end), 'DisplayName', 'Voltage and Hemo')
            hold on
            temp = lags(i,tempInd:end);
            xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

            legend('Location','northoutside')
            axis square
            xlabel('time (s)')
            ylabel('Correlation')
            title([ROIs(i+2).name(1:end-4)])
        end

        %save avg figure for each trial
        fname = ['Trial_Avg_', mouse,'_0',num2str(date),'_CrossCorr_Values'];
        saveFig(f4, fname, FolderName)
        fname = 'Scatter_Plot__Volt_Hemo_Temp';
        saveFig(f8, fname, FolderName)
        fname = 'Scatter_Plot__Volt_Pupil_Temp';
        saveFig(f9, fname, FolderName)
        fname = 'Scatter_Plot__Hemo_Pupil_Temp';
        saveFig(f10, fname, FolderName)
        fname = ['Trial_Avg_', mouse,'_0',num2str(date),'_CrossCorr_Values_VoltHemo'];
        saveFig(f11, fname, FolderName)
        close (f4, f11);

        %save animal specific data into larger Array
        allCorr = cat(3,allCorr, avgCorr);
        allLags = cat(3,allLags, avgLags);
        allCorrVP = cat(3,allCorrVP, avgCorrVP);
        allCorrHP = cat(3,allCorrHP, avgCorrHP);

    end

    % plot avg figure for all 5 animals 
    % Cross Corr
    f4 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:9
        subplot(2,5,i)
        tempCorr = mean(allCorr(i,:,:),3);
        [maxVal,maxInd] = max(abs(tempCorr));
        plot(lags(i,:),tempCorr, 'DisplayName', 'Voltage and Hemo')
        hold on
        axis square
        temp = lags(i,:);
        xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

        tempCorr = mean(allCorrVP(i,:,:),3);
        [maxVal,maxInd] = max(abs(tempCorr));
        plot(lags(i,:),tempCorr, 'DisplayName', 'Voltage and Pupil')
        hold on
        axis square
        temp = lags(i,:);
        xline(temp(maxInd),'r','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

        tempCorr = mean(allCorrHP(i,:,:),3);
        [maxVal,maxInd] = max(abs(tempCorr));
        plot(lags(i,:),tempCorr, 'DisplayName', 'Hemo and Pupil')
        hold on
        axis square
        temp = lags(i,:);
        xline(temp(maxInd),'y','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])
       
        legend('Location','northoutside')
        axis square
        xlabel('time (s)')
        ylabel('Correlation')
        title([ROIs(i+2).name(1:end-4)])
    end

    % Cross Corr Volt and Hemo in a small + window 
    f11 = figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:9
        subplot(2,5,i)
        tempInd = find(lags(i,:) == 0);
        startInd = tempInd - 250;
        endInd = tempInd + 250;
        tempCorr = mean(allCorr(i,:,:),3);
        [maxVal,maxInd] = min(tempCorr(tempInd:end));
        plot(lags(i,tempInd:end),tempCorr(tempInd:end), 'DisplayName', 'Voltage and Hemo')
        hold on
        axis square
        temp = lags(i,tempInd:end);
        xline(temp(maxInd),'b','DisplayName', ['Corr:', num2str(maxVal),' Delay:', num2str(temp(maxInd))])

        legend('Location','northoutside')
        axis square
        xlabel('time (s)')
        ylabel('Correlation')
        title([ROIs(i+2).name(1:end-4)])
    end

    %save avg figure for each trial
    fname = 'Avg_CrossCorr_Values';
    saveFig(f4, fname, FolderName)
    fname = 'Avg_CrossCorr_Values_VoltHemo';
    saveFig(f11, fname, FolderName)
    fname = 'Scatter_Plot__Volt_Hemo';
    saveFig(f8, fname, FolderName)
    fname = 'Scatter_Plot__Volt_Pupil';
    saveFig(f9, fname, FolderName)
    fname = 'Scatter_Plot__Hemo_Pupil';
    saveFig(f10, fname, FolderName)
    close all;
end

function [globalSignal] = calcGlobalSingal(data3D)
    [xDim, yDim, zDim] = size(data3D);
    data3D(isnan(data3D))=0;
    data2D = reshape(data3D, xDim * yDim, zDim);
    data2D = zscore(data2D, 0, 2);
    globalSignal = nanmean(data2D, 1);
end
function [gsRegressed] = GSR(signalData, globalSignal)
    
%     bft = signalData; % trace to be regressed: green channel
%     rg = globalSignal;  % trace to regress from: global signal
%     
%     % Calculate beta
%     beta=(rg'* bft)\(rg' * rg); 
% 
%     % regress  gs from data
%     gsRegressed = bft - rg*beta;

    %signalData = reshape(signalData, 10000, []);
    tmp = signalData'; 
    g = globalSignal';
    tmp = tmp-g*(g\tmp);
    gsRegressed = tmp';

end