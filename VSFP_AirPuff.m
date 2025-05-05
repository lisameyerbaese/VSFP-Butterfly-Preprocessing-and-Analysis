% Function loads in the preproccessed VSFP whisker stim data and plots the
% voltage and hemodynamic responses
%
%   Example Call VSFP_AirPuff('Cont', 'Right', 'Barrel')
%                VSFP_AirPuff('Sweep', 'Left', 'Barrel','Projection', 'yGSR')
%                VSFP_AirPuff('Sweep', 'Left', 'Barrel','Ratiometric', 'yGSR')
%
% Written by Lisa Meyer-Baese

function VSFP_AirPuff(Stim,Side, ROI, signalType,ynGSR)
    
    close all;
    %get root of file name based on computer being used 
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    folderName = [startFile,'\VSFP ButterFly\Data\VSFP_AirPuff\',Stim, '_',Side,'_', ROI,'_',signalType,'_',ynGSR];
    %check to see if Folder Exists to Save data, if not create it
     if ~exist(folderName, 'dir')
       mkdir(folderName)
    end

    allTrials = [dir([startFile,'\VSFP ButterFly\Data\VSFP_WhiskerStim\VSFP00_0917\', '*.mat']); dir([startFile,'\VSFP ButterFly\Data\VSFP_WhiskerStim\VSFP01_0918\', '*.mat'])];
    
    % load in BPOD Data for both mice
    fileBPOD = [];
    load([startFile, '\VSFP ButterFly\Data\VSFP_WhiskerStim\BPOD\','VSFP00_0917_BPOD.mat']) % has data for all 100 trials
    fileBPOD = [fileBPOD, allFiles];
    load([startFile, '\VSFP ButterFly\Data\VSFP_WhiskerStim\BPOD\','VSFP01_0918_BPOD.mat']) % has data for only the first 50 trials
    fileBPOD = [fileBPOD, allFiles];

    % load in whole cortical map
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));
    fixed_map = imresize(fixed_map, [100, 100]);
    M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts

    fs = 200;
    T = 1/fs;             % Sampling period       
    L = 4096;             % Length of signal
    time = (0:L-1)*T;     % Time vector

    roiType = 0;
    
    allVoltOG = [];
    allHemoOG = [];
    allDonor = [];
    allAcceptor = [];
    allGSCrossCorr = [];
    allDiffGS = [];
    

    for i = 1:length(allTrials)
        
        %load in information about the stimulus type and stimulus side
        %load([allTrials(i).folder,'\',allTrials(i).name], 'stimType', 'stimSide')
        clear stimSide
        load([allTrials(i).folder,'\',allTrials(i).name],'stimSide')
        if ~exist('stimSide','var')
            load([allTrials(i).folder,'\',allTrials(i).name])
            if str2double(out.fNum) < 50 
                stimSide = 'Left';
            else
                stimSide = 'Right';
            end
            save([allTrials(i).folder,'\',allTrials(i).name], 'out','stimSide')
        end
        
        %if isequal(stimType,Stim) && isequal(stimSide, Side)
        if isequal(stimSide, Side)
            
            % load in the preProcessed imaging data using Rhett's code
            load([allTrials(i).folder,'\',allTrials(i).name])
            mouseID = out.mouseID;
            trial = out.fNum;
            fileName = out.fileNum;

            % get the correct start of imaging for that trial
            tempTrials = [fileBPOD.trialNum];
            index = ismember({fileBPOD.mouseID},mouseID) & ismember(tempTrials,str2double(trial)) & ismember({fileBPOD.StimSide},stimSide) ;
            delaySec = fileBPOD(index).ImgSyncTriggerTimes - 2;
            delayFrames = round(delaySec * fs);

            % use avg green channel to select Barrel ROI
            if roiType == 0 && isequal(mouseID,'VSFP00')

                % see if saved ROI already exists
                if isfile([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_','_ROI.mat'])
%                     load([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_','_ROI.mat'])
%                 else
                    f1 = figure(1);
                    meanImg = mean(out.imgD,3);
                    axis square
                    I = mat2gray(meanImg);
                    imshow(I, 'InitialMagnification','fit');
                    title([mouseID, ' Barrel ROI ', Stim,' ',Side])
                    r1 = drawcircle('Label','ROI','Color',[1 0 0]);
                    maskImg = createMask(r1);
                    roiType = 1;
                    save([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_','_ROI.mat'],'maskImg')
                    fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
                    saveFig(f1, [fname, '_ROI'],folderName)
                    close(f1)
                end

            elseif roiType == 1 && isequal(mouseID,'VSFP01')
                
                 % see if saved ROI already exists
                if isfile([allTrials(i).folder,'\', mouseID,'_ROI.mat'])
%                     load([allTrials(i).folder,'\', mouseID,'_ROI.mat'])
%                 else
                    roiType = 0;
                    f1 = figure(1);
                    meanImg = mean(out.imgD,3);
                    I = mat2gray(meanImg);
                    imshow(I, 'InitialMagnification','fit');
                    r1 = drawcircle('Label','ROI','Color',[1 0 0]);
                    title([mouseID, ' Barrel ROI ', Stim,' ',Side])
                    center1 = r1. Center; 
                    radius1 = r1.Radius;
                    maskImg = createMask(r1);
                    save([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_','_ROI.mat'],'maskImg')
                    fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
                    saveFig(f1, [fname, '_ROI'],folderName)
                    close(f1)
                end
            end
            
            % make a plot of the scatter plot of raw donor and acceptor,
            % baseline is first 200ms
            rawDonor = out.imgD;
            rawAcceptor = out.imgA;

            [xDim, yDim, zDim] = size(rawDonor);
            avgingFr = floor(fs./5);
            donor2D = reshape(rawDonor, xDim * yDim, zDim);
            baselineDonor = mean(donor2D(:,1:avgingFr),2);
            dfDonor = (donor2D - baselineDonor) ./ baselineDonor;
            dfDonor1D = reshape(dfDonor, 1, []);
        
            acceptor2D = reshape(rawAcceptor, xDim * yDim, zDim);
            baselineAcceptor = mean(acceptor2D(:,1:avgingFr),2);
            dfAcceptor = (acceptor2D - baselineAcceptor)./ baselineAcceptor;
            dfAcceptor1D = reshape(dfAcceptor, 1, []);
            
            randInd = randsample(1:40960000,1000000);
           
            f1 = figure(1);
            s = scatter(dfAcceptor1D(randInd),dfDonor1D(randInd) ,'k','.');
            s.AlphaData = 0.5;
            axis equal
            axis square
            xlabel('Acceptor dF/F (%)')
            ylabel('Donor dF/F (%)')
            title('Scatter Plot of Raw Data')
            fname = [allTrials(i).name(1:end-4)];
            saveFig(f1, [fname,'_DonorAcceptor_ScatterPlot'],folderName)
            close(f1)

            % compare Donor/Acceptor Data to voltage and hemodynamic
            % response
            f2 = figure(2);
            subplot(3,2,1)
            dROI = out.imgD .* maskImg;
            aROI = out.imgA .* maskImg;
            dROI2D = reshape(dROI, 10000, []);
            aROI2D = reshape(aROI, 10000, []);
            plot(time, zscore(mean(dROI2D,1)),'g', 'DisplayName', 'Donor')
            hold on 
            plot(time, zscore(mean(aROI2D,1)), 'r','DisplayName', 'Acceptor')
            title('Donar/Acceptor Response')
            legend()
            hold on
            addLines(delaySec)

            subplot(3,2,2)
            temp = zscore(mean(dROI2D,1));
            plot(time(801:1505),temp(801:1505) ,'g', 'DisplayName', 'Donor')
            hold on 
            temp = zscore(mean(aROI2D,1));
            plot(time(801:1505),temp(801:1505) , 'r','DisplayName', 'Acceptor')
            title('Donar/Acceptor Response')
            xline(5 - delaySec)
            xline(6 - delaySec)
            xline(7 - delaySec)
            legend()
            clear temp

            % load in volt and hemo 
            if isequal(signalType,'Projection')
                maskVolt = out.projectedVolt .* maskImg;
                maskVolt2D = reshape(maskVolt, 10000,[]);
                maskVolt2D(isnan(maskVolt2D))=0;
                maskHemo = out.projectedHemo .* maskImg;
                maskHemo2D = reshape(maskHemo, 10000, []);
                maskHemo2D(isnan(maskHemo2D))=0;
                if isequal(ynGSR, 'yGSR')
                    lpHemo2D = lowpass(maskHemo2D',5,vsfp_data.out.Fs); %filters each column independantly
                    hemo = reshape(lpHemo2D',100,100,[]);
                    tempMask = imrotate(vsfp_data.out.mask, 270);
                    gsHemo = calcGlobalSingal(hemo.* tempMask);
                    [maskHemo2D] = GSR(maskHemo2D, gsHemo);
                end
              
            else
                maskVolt = out.imgDR3 .* maskImg;
                maskVolt2D = reshape(maskVolt, 10000,[]);
                maskHemo = out.imgDS .* maskImg;
                maskHemo2D = reshape(maskHemo, 10000, []);
                gsHemo = out.gsHemo;
                if isequal(ynGSR, 'yGSR')
                    [maskHemo2D] = GSR(maskHemo2D, gsHemo);
                end
            end
            
            subplot(3,2,3)
            
            plot(time, zscore(mean(maskVolt2D,1)),'k', 'DisplayName', 'Volt')
            hold on 
            title('Voltage Response')
            addLines(delaySec)

            subplot(3,2,4)
            temp = zscore(mean(maskVolt2D,1));
            plot(time(801:1505),temp(801:1505) ,'k', 'DisplayName', 'Volt')
            hold on
            xline(5 - delaySec)
            xline(6 - delaySec)
            xline(7 - delaySec)
            title('Volt')

            subplot(3,2,5)
            plot(time, zscore(mean(maskHemo2D,1)),'g', 'DisplayName', 'Hemo')
            hold on 
            title('Hemodynamics Response')
            hold on
            addLines(delaySec)

            subplot(3,2,6)
            temp = zscore(mean(maskHemo2D,1));
            plot(time(801:1505),temp(801:1505) ,'g', 'DisplayName', 'Volt')
            hold on 
            xline(5 - delaySec)
            xline(6 - delaySec)
            xline(7 - delaySec)
            title('Hemo')

            % look at correlation between global singal for volt/hemo on raw
            % and on run smoothed traces by 1 sec
            f3 = figure(3);
            tempHemo = mean(maskHemo2D,1);
            smoothHemo = movmean(tempHemo,fs * 1);
            subplot(3,1,1)
            plot(time, tempHemo, 'g', 'DisplayName', 'Hemo')
            addLines(delaySec)
            hold on 
            plot(time, smoothHemo, 'color',"#77AC30", 'DisplayName', 'Smooth Hemo')
            tempVolt = mean(maskVolt2D,1);
            smoothVolt = movmean(tempVolt, fs * 1);
            plot(time, tempVolt, 'color',[.7 .7 .7], 'DisplayName', 'Volt')
            plot(time, smoothVolt,'color', 'k', 'DisplayName', 'Smooth Volt')
            xlabel('Time (s)')
            legend()

            % cross correlation between original traces
            subplot(3,1,2)
            [cOG,lags] = xcorr(tempHemo,tempVolt,'normalized');
            lags = lags / fs;
            plot(lags,cOG)
            xlabel('Time (s)')
            title('Cross Corr OG Hemo w/ Volt')
            
            % correlation between smoothed traces
            subplot(3,1, 3)
            [cOG,lags] = xcorr(smoothHemo,smoothVolt,'normalized');
            lags = lags / fs;
            plot(lags,cOG)
            xlabel('Time (s)')
            title('Cross Corr Smoothed Hemo w/ Volt')
           
            fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
            saveFig(f2, [fname, '_RawSignals'],folderName)
            saveFig(f3, [fname, '_CrossCorr'],folderName)
            close(f2,f3)

            % look at the cross corr between GS hemo and GS Volt
            f10 = figure(10);
            gsVolt = out.gsVolt;
            gsHemo = -1 * out.gsHemoLP;
            subplot(3,1,1)
            plot(time, gsVolt, 'k')
            hold on
            [~,indVolt] = max(gsVolt);
            xline(time(indVolt))
            title(['Global Volt max at t=', num2str(time(indVolt))])
            subplot(3,1,2)
            plot(time, gsHemo, 'g')
            [~,indHemo] = max(gsHemo);
            xline(time(indHemo))
            title(['Inverted Global Hemo max at t=', num2str(time(indHemo))])
            Fs = out.Fs;             % Sampling Frequency (Hz)

            % Filter the hemo signal at 0.1Hz prior to calculating the
            % cross corr
            gsHemoLP = lowpass(gsHemo,1,Fs,ImpulseResponse="iir",Steepness=0.95);
            [r,lags] = xcorr(zscore(gsHemoLP),zscore(gsVolt),'normalized', Fs*5);
            t = lags/Fs;             % Time Vector (seconds)
            subplot(3,1,3)
            plot(t, r)
            title('Cross Corr (hemo/volt)')
            
            fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
            saveFig(f10, [fname, '_GS_CrossCorr'],folderName)
            close(f10)

            f11 = figure(11);
            grayColor = [.7 .7 .7];
            plot(t, r,'--', 'Color', grayColor)
            hold on
            title('Cross Corr Global (hemo/volt)')

            % save global signal cross corr values and the difference
            % between max peaks 
            allGSCrossCorr = [allGSCrossCorr; r];
            allDiffGS = [allDiffGS; (time(indHemo) - time(indVolt))];

            % save traces and correlations to look at group avg based on
            % individual delay, aligned traces to onset of first airPuff
            % shows 2 sec of baseline + all airpuff traces (15s) + 0.5 sec
            % baseline
            indTemp = (1000-delayFrames) - 400:1:(1000-delayFrames) + 3100;
            temp = zscore(mean(dROI2D,1));
            allDonor = [allDonor;temp(indTemp)];
            temp = zscore(mean(aROI2D,1));
            allAcceptor = [allAcceptor;temp(indTemp)];
            temp = zscore(mean(maskVolt2D,1));
            allVoltOG = [allVoltOG;temp(indTemp)];
            temp = zscore(mean(maskHemo2D,1));
            allHemoOG = [allHemoOG; temp(indTemp)];
            clear indTemp

        end
    end
    
    % plot the average activity as a funciton of time
    f5 = figure(5);
    timeTemp = [-(Fs * 2):1: (15.5 * Fs)]/Fs;
    subplot(3,2,1)
    tempVolt = mean(allDonor,1);
    plot(timeTemp, tempVolt, 'g')
    hold on 
    tempHemo = mean(allAcceptor,1);
    plot(timeTemp, tempHemo, 'r')
    addLines(5)
    xlabel('Time (s)')
    title('Avg Donor/Acceptor Traces')

    subplot(3,2,2)
    timeTempShort = [-(Fs * 2):1: (3 * Fs)-1]/Fs;
    plot(timeTempShort, tempVolt(1:1:5*Fs), 'g')
    hold on 
    plot(timeTempShort, tempHemo(1:1:5*Fs), 'r')
    xlabel('Time (s)')
    title('Zoomed in Avg Donor/Acceptor Traces')
    
    subplot(3,2,3)
    tempVolt = mean(allVoltOG,1);
    plot(timeTemp,tempVolt, 'k')
    hold on 
    tempHemo = mean(allHemoOG,1);
    plot(timeTemp,tempHemo, 'g')
    addLines(5)
    xlabel('Time (s)')
    title('Avg Traces')

    subplot(3,2,4)
    tempVolt = mean(allVoltOG,1);
    plot(timeTempShort,tempVolt(1:1:5*Fs), 'k')
    hold on 
    tempHemo = mean(allHemoOG,1);
    plot(timeTempShort,tempHemo(1:1:5*Fs), 'g')
    %addLines(3)
    xlabel('Time (s)')
    title('Avg Traces')
    
    subplot(3,2,5)
    [cTemp,lTemp] = xcorr(tempHemo,tempVolt, fs* 10,'normalized');
    lTemp = lTemp /fs;
    plot(lTemp,cTemp)
    xlabel('Time (s)')
    title('Avg Cross Corr ')

    % plot the fft of the voltage and hemo activity for different Freq Stim
    Fs = 200;            % Sampling frequency                    
    T = 1/Fs;            % Sampling period       
    L = 3 * Fs;          % Length of signal

    xVolt = mean(allVoltOG,1);
    xHemo = mean(allHemoOG,1);
    tPointsStart = [2,8,14];
    tPointsEnd = [8,14,17];
    stimFreq = [1,2,5];
    
    % look at activity for each of the stim frequencies
    for freq = 1:length(tPointsStart)
        startT = tPointsStart(freq);
        endT = tPointsEnd(freq);

        tempVolt = xVolt(1,startT*200:endT*200);
        tempHemo = xHemo(1,startT*200:endT*200);

        f6 = figure('units','normalized','outerposition',[0 0 1 1]);
        % fft Volt
        yVolt = fft(tempVolt(1:end-1));   
        P2 = abs(yVolt/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs/L*(0:(L/2));
        %truncate the x-axis to start at 0.667Hz
        subplot(2,1,1)
        plot(f(1,3:end),P1(1,3:end),"LineWidth",3,'Color','k')
        %plot(f,P1,"LineWidth",3,'Color','k')
        maxVolt = max(P1(1,3:end));
        hold on
        xline(stimFreq(freq))
        title('Avg Volt in Barrel')
        xlabel("f (Hz)")
        ylabel("Amplitude")
    
        % fft hemo
        yHemo = fft(tempHemo(1:end-1));
        P2 = abs(yHemo/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs/L*(0:(L/2));

        %truncate the x-axis to start at 0.667Hz
        subplot(2,1,2)
        plot(f(1,3:end),P1(1,3:end),"LineWidth",3, 'Color',"#77AC30") 
        %plot(f,P1,"LineWidth",3, 'Color',"#77AC30")
        maxHemo = max(P1(1,3:end));
        hold on
        xline(stimFreq(freq))
        title('Avg Hemo in Barrel')
        xlabel("f (Hz)")
        ylabel("Amplitude")
        sgtitle([num2str(stimFreq(freq)),' Hz']);
        %ylim([0 max(maxHemo, maxVolt)])

        subplot(2,1,1)
        ylim([0 max(maxHemo, maxVolt)])
        saveFig(f6, ['Avg_',num2str(stimFreq(freq)),'_',Stim,'_', Side, '_VoltHemo_FFT'],folderName)
        close(f6)
        
    end

    saveFig(f5, ['Avg_',Stim,'_', Side, '_VoltHemo_CrossCorr'],folderName)
    save([folderName, '\AvgActivity.mat'],'allDonor','allAcceptor','allVoltOG','allHemoOG')

    % plot the average GS correlation 
    f7 = figure(7);
    plotAvgCrossCorr(lags, Fs, allGSCrossCorr, 1, 1, 'gsHemo', 'gsVolt')
    sgtitle(['Avg Delay', num2str(mean(allDiffGS))])
    saveFig(f7, ['Avg_',Stim,'_', Side, '_GS_CrossCorr'],folderName)

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

function addLines(delaySec)
    % this function adds the lines to indicate the start and end of stimulation
    
    % for 1 Hz
    xline(5 - delaySec, 'Color', '#A2142F','HandleVisibility','off')
    xline(6 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(7 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')

    xline(8 - delaySec, 'Color', '#A2142F','HandleVisibility','off')
    xline(9 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(10 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')

    % for 2 Hz
    xline(11 - delaySec, 'Color', '#A2142F','HandleVisibility','off')
    xline(11.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(12 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(12.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(13 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(13.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')

    xline(14 - delaySec, 'Color', '#A2142F', 'HandleVisibility','off')
    xline(14.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(15 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(15.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(16 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(16.5 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')

    % for 5 Hz
    xline(17 - delaySec, 'Color', '#A2142F', 'HandleVisibility','off')
    xline(17.2 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(17.4 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(17.6 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(17.8 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(18 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(18.2 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(18.4 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(18.6 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(18.8 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(19 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(19.2 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(19.4 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(19.6 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(19.8 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')
    xline(20 - delaySec, 'Color', '#4DBEEE','HandleVisibility','off')

    
end
function [globalSignal] = calcGlobalSingal(data3D)
    [xDim, yDim, zDim] = size(data3D);
    data3D(isnan(data3D))=0;
    data2D = reshape(data3D, xDim * yDim, zDim);
    data2D = zscore(data2D, 0, 2);
    globalSignal = nanmean(data2D, 1);
end

function [gsRegressed] = GSR(signalData, globalSignal)

    tmp = signalData'; 
    g = globalSignal';
    tmp = tmp-g*(g\tmp);
    gsRegressed = tmp';

end