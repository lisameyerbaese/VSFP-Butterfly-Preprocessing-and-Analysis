% Function loads in the preproccessed VSFP whisker stim data and plots the
% voltage and hemodynamic responses and looks to see across different ROIs
% if the increase in light intensity which is time locked to stimulus onset is local
% or global
%
%   Example Call 
%                VSFP_AirPuff_TestMotion('Sweep', 'Left', 'Several','Projection')
%                VSFP_AirPuff_TestMotion('Sweep', 'Left', 'Several','Ratiometric')
%
% Written by Lisa Meyer-Baese

function VSFP_AirPuff_TestMotion(Stim,Side, ROI, signalType)
    
    close all;
    %get root of file name based on computer being used 
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    folderName = [startFile,'\VSFP ButterFly\Data\VSFP_AirPuff_TestMotion\',Stim, '_',Side,'_', ROI,'_',signalType];
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
 
    tempLen = length(allTrials);
    allVoltOG = cell(6,tempLen);
    allHemoOG = cell(6,tempLen);
    allDonor = cell(6,tempLen);
    allAcceptor = cell(6,tempLen);
    allgsVolt = cell(1,tempLen);
    allgsHemo = cell(1,tempLen);
    
    allVolt3D = [];
    allHemo3D = [];

    for i = 1:length(allTrials)
        
        %load in information about the stimulus type and stimulus side
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

            % use avg green channel to select 5 ROIs
            if roiType == 0 && isequal(mouseID,'VSFP00')
                f1 = figure(1);
                meanImg = mean(out.imgD,3);
                axis square
                I = mat2gray(meanImg);
                imshow(I, 'InitialMagnification','fit');
                title([mouseID, ' ROIs: Contra Barrel, Ipsi Barrel, Contra V1, Contra M1, and Contra Outside for ', Stim,' ',Side])
                r1 = drawcircle('Label','Contra Barrel','Color',[1 0 0]);
                mask1 = createMask(r1);
                r2 = drawcircle('Label','Ipsi Barrel','Color',[1 0 0]);
                mask2 = createMask(r2);
                r3 = drawcircle('Label','V1','Color',[1 0 0]);
                mask3 = createMask(r3);
                r4 = drawcircle('Label','M1','Color',[1 0 0]);
                mask4 = createMask(r4);
                r5 = drawcircle('Label','Outside','Color',[1 0 0]);
                mask5 = createMask(r5);
                allMasks = {mask1, mask2, mask3, mask4, mask5};
                roiType = 1;
                save([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_',ROI,'_ROI.mat'],'allMasks')
                fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
                saveFig(f1, [fname, '_ROI'],folderName)
            elseif roiType == 1 && isequal(mouseID,'VSFP01')
                roiType = 0;
                f1 = figure(1);
                meanImg = mean(out.imgD,3);
                I = mat2gray(meanImg);
                imshow(I, 'InitialMagnification','fit');
                title([mouseID, ' ROIs: Contra Barrel, Ipsi Barrel, Contra V1, Contra M2, and Contra Outside for ', Stim,' ',Side])
                r1 = drawcircle('Label','Contra Barrel','Color',[1 0 0]);
                mask1 = createMask(r1);
                r2 = drawcircle('Label','Ipsi Barrel','Color',[1 0 0]);
                mask2 = createMask(r2);
                r3 = drawcircle('Label','V1','Color',[1 0 0]);
                mask3 = createMask(r3);
                r4 = drawcircle('Label','M2','Color',[1 0 0]);
                mask4 = createMask(r4);
                r5 = drawcircle('Label','Outside','Color',[1 0 0]);
                mask5 = createMask(r5);
                allMasks = {mask1, mask2, mask3, mask4, mask5};
                save([allTrials(i).folder,'\', mouseID,'_',Stim,'_',Side,'_',ROI,'_ROI.mat'],'allMasks')
                fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
                saveFig(f1, [fname, '_ROI'],folderName)
            end
            

            % loop through all 5 ROIs and plot the donor/acceptor and
            % volt/hemo data
            f2 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(6,2,1)
            imshow(I, 'InitialMagnification','fit');
            hold on 
            viscircles(r1.Center, r1.Radius, 'Color',"#D95319")
            viscircles(r2.Center, r2.Radius, 'color',"#EDB120")
            viscircles(r3.Center, r3.Radius, 'color',"#0072BD")
            viscircles(r4.Center, r4.Radius, 'Color',"#77AC30")
            viscircles(r5.Center, r5.Radius, 'Color',"#4DBEEE")
            roiColors = {"#D95319","#EDB120","#0072BD","#77AC30","#4DBEEE"};
            roiLables = {'Contra Barrel', 'Ipsi Barrel', 'V1', 'M2', 'Outside'};
            axis off
            tempInd = 3;
            
            for r = 1:5
                maskImg = allMasks{r};
                % get raw data from donor and acceptor
                dROI = out.imgD .* maskImg;
                aROI = out.imgA .* maskImg;
                dROI2D = reshape(dROI, 10000, []);
                aROI2D = reshape(aROI, 10000, []);
    
                % load in volt and hemo 
                if isequal(signalType,'Projection')
                    maskVolt = out.projectedVolt .* maskImg;
                    maskVolt2D = reshape(maskVolt, 10000,[]);
                    maskVolt2D(isnan(maskVolt2D))=0;
                    maskHemo = out.projectedHemo .* maskImg;
                    maskHemo2D = reshape(maskHemo, 10000, []);
                    maskHemo2D(isnan(maskHemo2D))=0; 
                    volt3D = out.projectedVolt;
                    hemo3D = out.projectedHemo;
                else
                    maskVolt = out.imgDR3 .* maskImg;
                    maskVolt2D = reshape(maskVolt, 10000,[]);
                    maskHemo = out.imgDS .* maskImg;
                    maskHemo2D = reshape(maskHemo, 10000, []);
                    volt3D = out.imgDR3;
                    hemo3D = out.imgDS;
                end
                
                % compare Donor/Acceptor Data to voltage and hemodynamic
                % response for different ROIs
                subplot(6,2,tempInd)
                temp = zscore(mean(dROI2D,1));
                plot(time(801:1505),temp(801:1505) ,'g')
                hold on 
                temp = zscore(mean(aROI2D,1));
                plot(time(801:1505),temp(801:1505) , 'r')
                title('Donor (green) and Acceptor (red) Response', 'Color', roiColors{r})
                xline(5 - delaySec)
                xline(6 - delaySec)
                xline(7 - delaySec)
                ylabel(['Z-score ' roiLables{r}])
                xlabel('Time(s)')
                clear temp

                subplot(6,2,tempInd+1)
                temp = zscore(mean(maskVolt2D,1));
                plot(time(801:1505),temp(801:1505) ,'k')
                hold on
                xline(5 - delaySec)
                xline(6 - delaySec)
                xline(7 - delaySec)
                hold on
                temp = zscore(mean(maskHemo2D,1));
                plot(time(801:1505),temp(801:1505) ,'color',"#0072BD")
                hold on 
                xline(5 - delaySec)
                xline(6 - delaySec)
                xline(7 - delaySec)
                title('Volt (black) and Hemo (blue)', 'Color', roiColors{r})
                ylabel(['Z-score ' roiLables{r}])
                xlabel('Time(s)')

                % save individual traces for Trial Avg, 2 sec baseline plus
                % all the stim data
                indTemp = (1000-delayFrames) - 400:1:(1000-delayFrames) + 3100;
                temp = zscore(mean(dROI2D,1));
                allDonor{r,i} = temp(indTemp);
                temp = zscore(mean(aROI2D,1));
                allAcceptor{r,i} = temp(indTemp);
                temp = zscore(mean(maskVolt2D,1));
                allVoltOG{r,i} = temp(indTemp);
                temp = zscore(mean(maskHemo2D,1));
                allHemoOG{r,i} = temp(indTemp);
                tempInd = tempInd + 2;

            end
            
            subplot(6,2,2)
            % add the GS for this trial to the plot
            gsVolt = out.gsVolt(801:1505);
            gsHemo = out.gsHemoLP(801:1505);
            plot(time(801:1505), gsVolt, 'k')
            hold on
            plot(time(801:1505), gsHemo,'color',"#0072BD")
            title('Global Volt (black) and Hemo (blue)')
            hold on
            ylabel('Z-score GS')
            xline(5 - delaySec)
            xline(6 - delaySec)
            xline(7 - delaySec)
            figTitle = regexprep( [allTrials(i).name(1:end-4),'_',Stim,'_',Side], '_', ' ');
            sgtitle(['Activity Across ROIs for ', figTitle])
            
            % make figure with maps of imaging data for specific ranges
            f3 = figure('units','normalized','outerposition',[0 0 1 1]);
            % plot imaging data for 5 periods in time: avg of baseline, 0-0.02 post
            % stim, 0.085-0.105 post stim, and 0.5, 1 and 2s post stim. Each plot is
            % avg of 0.02 sec, or 4 samples
            range1 = 410:1:414;
            range2 = 417:1:421;
            range3 = 500:1:504;
            range4 = 600:1:604;
            range5 = 800:1:804;
            
            temp = volt3D(:,:,indTemp) .* out.mask;
            baseVal = mean(temp(:,:,1:1:200),3);

            subplot(2,6,1)
            imagesc(imrotate(baseVal, 270), [-0.1, 0.1]), colorbar
            axis square, axis off
            title('Baseline')
        
            subplot(2,6,2)
            stimVal = (mean(temp(:,:,range1),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.05 - 0.07s post stim')
        
            subplot(2,6,3)
            stimVal = (mean(temp(:,:,range2),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.085 - 0.1050s post stim')
        
            subplot(2,6,4)
            stimVal = (mean(temp(:,:,range3),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.5 - 0.052s post stim')
        
            subplot(2,6,5)
            stimVal = (mean(temp(:,:,range4),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('1 - 1.02s post stim')
        
            subplot(2,6,6)
            stimVal = (mean(temp(:,:,range5),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('2 - 2.02s post stim')
            
            temp = hemo3D(:,:,indTemp) .* out.mask;
            baseVal = mean(temp(:,:,1:1:400),3);
           
            subplot(2,6,7)
            imagesc(imrotate(baseVal, 270)), colorbar
            axis square, axis off
            title('Baseline')
        
            subplot(2,6,8)
            stimVal = (mean(temp(:,:,range1),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.05 - 0.07s post stim')
        
            subplot(2,6,9)
            stimVal = (mean(temp(:,:,range2),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.085 - 0.1050 s post stim')
        
            subplot(2,6,10)
            stimVal = (mean(temp(:,:,range3),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('0.5 - 0.052 s post stim')
        
            subplot(2,6,11)
            stimVal = (mean(temp(:,:,range4),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('1 - 1.02s post stim')
            
            subplot(2,6,12)
            stimVal = (mean(temp(:,:,range5),3));
            stimVal = imrotate(stimVal, 270);
            imagesc(stimVal), colorbar
            axis square, axis off
            title('2 - 2.02s post stim')
            sgtitle([figTitle ' Top:Voltage Bottom: OG Hemo (not sign inverted). Not Relative to Baseline'])
            
            % save GS signals 
            allgsVolt{1,i} = zscore(out.gsVolt(indTemp));
            allgsHemo{1,i} = zscore(out.gsHemoLP(indTemp));

            % save whole frames as well to spatial maps
            allVolt3D = cat(4,allVolt3D, volt3D(:,:,indTemp) .* out.mask);
            allHemo3D = cat(4, allHemo3D, hemo3D(:,:,indTemp) .* out.mask);

            fname = [allTrials(i).name(1:end-4),'_',Stim,'_',Side];
            saveFig(f2, [fname, '_ROI_Traces'],folderName)
            saveFig(f3, [fname, '_Spatial_Maps'],folderName)
            close(f2, f3)
        end
    end
    
    % plot the average activity as a funciton of time
    f5 = figure('units','normalized','outerposition',[0 0 1 1]);
%     timeTemp = [-(fs * 2):1:(3 * fs)]/fs;
%     indTemp = 1:1:1001;
    timeTemp = [-(fs * 2):1:(15.505 * fs)]/fs;
    timeTemp(end) = [];
    indTemp = 1:1:3501;
    subplot(6,2,1)
    imshow(I, 'InitialMagnification','fit');
    hold on 
    viscircles(r1.Center, r1.Radius, 'Color',"#D95319")
    viscircles(r2.Center, r2.Radius, 'color',"#EDB120")
    viscircles(r3.Center, r3.Radius, 'color',"#0072BD")
    viscircles(r4.Center, r4.Radius, 'Color',"#77AC30")
    viscircles(r5.Center, r5.Radius, 'Color',"#4DBEEE")
    roiColors = {"#D95319","#EDB120","#0072BD","#77AC30","#4DBEEE"};
    roiLables = {'Contra Barrel', 'Ipsi Barrel', 'V1', 'M1', 'Outside'};
    axis off

    % need to add xlines at right spots and plot the avg for other ROIs
    subplot(6,2,2)
    allGSVolt = mean(cell2mat(allgsVolt'),1);
    allGSHemo = mean(cell2mat(allgsHemo'),1);
    plot(timeTemp, allGSVolt(indTemp), 'k')
    hold on 
    plot(timeTemp, allGSHemo(indTemp), 'b')
    xlabel('Time (s)')
    ylabel('z-scored')
    title('Avg GS Traces')
    xline(0); xline(1); xline(2); xline(3); xline(4); xline(5);
    ylim([-1,1])
    %axis square

    % plot donor/acceptor and volt/hemo avg for all other ROIs
    indPlot = 3;
    for r = 1:5
        subplot(6,2,indPlot)
        tempDonor = mean(cell2mat(allDonor(r,:)'),1);
        tempAcceptor = mean(cell2mat(allAcceptor(r,:)'),1);
        plot(timeTemp,tempDonor(indTemp),'g')
        hold on 
        plot(timeTemp,tempAcceptor(indTemp),'r')
        title('Donar (green) and Acceptor (red) Response', 'Color', roiColors{r})
        ylabel('z-score')
        xlabel('Time(s)')
        xline(0); xline(1); xline(2); xline(3); xline(4); xline(5);
        ylim([-1,1])
        %axis square

        subplot(6,2,indPlot+1)
        tempVolt = mean(cell2mat(allVoltOG(r,:)'),1);
        tempHemo = mean(cell2mat(allHemoOG(r,:)'),1);
        plot(timeTemp,tempVolt(indTemp),'k')
        hold on 
        plot(timeTemp,tempHemo(indTemp),'r','color',"#0072BD")
        title('Volt (black) and Hemo (blue)', 'Color', roiColors{r})
        ylabel('z-score')
        xlabel('Time(s)')
        xline(0); xline(1); xline(2); xline(3); xline(4); xline(5);
        ylim([-1,1])
        %axis square
        
        indPlot = indPlot +2;
    end
    sgtitle(['Trial Avg for ' Stim ' ' Side])
    saveFig(f5, ['Avg_',Stim,'_', Side, '_Avg_Traces'],folderName)
    close(f5)

    % plot imaging data for 5 periods in time: avg of baseline, 0-0.02 post
    % stim, 0.085-0.105 post stim, and 0.5, 1 and 2s post stim. Each plot is
    % avg of 0.02 sec, or 4 samples
    f6 = figure(6);
    range1 = 410:1:414;
    range2 = 417:1:421;
    range3 = 500:1:504;
    range4 = 600:1:604;
    range5 = 800:1:804;
    
    temp = mean(allVolt3D,4);
    baseVal = mean(temp(:,:,1:1:400),3);

    subplot(2,6,1)
    imagesc(imrotate(baseVal, 270) * 100, [-1 1])
    axis square, axis off
    title('Baseline')

    subplot(2,6,2)
    stimVal = (mean(temp(:,:,range1),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal , [-1 1])
    axis square, axis off
    title('0.05 - 0.07s post stim')

    subplot(2,6,3)
    stimVal = (mean(temp(:,:,range2),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal , [-1 1])
    axis square, axis off
    title('0.085 - 0.1050s post stim')

    subplot(2,6,4)
    stimVal = (mean(temp(:,:,range3),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal , [-1 1])
    axis square, axis off
    title('0.5 - 0.052s post stim')

    subplot(2,6,5)
    stimVal = (mean(temp(:,:,range4),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal , [-1 1])
    axis square, axis off
    title('1 - 1.02s post stim')

    subplot(2,6,6)
    stimVal = (mean(temp(:,:,range5),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('2 - 2.02s post stim')
    
    temp = mean(allHemo3D,4);
    baseVal = mean(temp(:,:,1:1:400),3);
   
    subplot(2,6,7)
    imagesc(imrotate(baseVal, 270)  * 100, [-1 1])
    axis square, axis off
    title('Baseline')

    subplot(2,6,8)
    stimVal = (mean(temp(:,:,range1),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('0.05 - 0.07s post stim')

    subplot(2,6,9)
    stimVal = (mean(temp(:,:,range2),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('0.085 - 0.1050 s post stim')

    subplot(2,6,10)
    stimVal = (mean(temp(:,:,range3),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('0.5 - 0.052 s post stim')

    subplot(2,6,11)
    stimVal = (mean(temp(:,:,range4),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('1 - 1.02s post stim')
    
    subplot(2,6,12)
    stimVal = (mean(temp(:,:,range5),3) - baseVal) * 100;
    stimVal = imrotate(stimVal, 270);
    imagesc(stimVal, [-1 1])
    axis square, axis off
    title('2 - 2.02s post stim')
    sgtitle('Trial Avg Response Top:Voltage while Bottom: OG Hemo (not sign inverted). All are baseline subtracted')

    save([folderName, '\AvgActivity.mat'],'allDonor','allAcceptor','allVoltOG','allHemoOG', 'allGSHemo','allGSVolt')
end


