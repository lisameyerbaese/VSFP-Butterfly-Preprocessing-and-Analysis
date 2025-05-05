% Function takes xlsx sheet and generates figures for following analysis
%   1) Takes raw traces of pupil,volt and hemo signal
%   Creates Magnitude Scalogram and Wavelet Coherence Plots for Pairs of
%   those signals
%
% Inputs: 
%        VSFP_WaveletCoherence_ROI(ROI) to specifiy which region to look at
%        specifically 
%
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%
% Written by Lisa Meyer-Baese
%
% update 11/8/23 print coherence in semi-log scale to make low freq easier
% to resolve

function VSFP_WaveletCoherence_ROI(region)
    close all;
    
    % file directory and color properties depending on which ROI is choosen
    switch region
        case 1
            ROI = 'VSFP_WaveletROI_one\';
            color = load('YIOrBr8.mat');
            color = color.YlOrBr8;
        case 2
            ROI = 'VSFP_WaveletROI_two\';
            color = load('Greens8.mat');
            color = color.Greens8;
        case 3
            ROI = 'VSFP_WaveletROI_three\';
            color = load('Blues8.mat');
            color = color.Blues8;
        case 4
    
            ROI = 'VSFP_WaveletROI_four\';
            color = load('Greys8.mat');
            color = color.Greys8;
    end
            

    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath 'C:\Users\lmeyerb\OneDrive - Emory University\Documents\MATLAB\Examples\R2023b\wavelet\WaveletCoherenceAnalysisOfBrainDynamicsExample'
    addpath([startFile, '\VSFP ButterFly\Info'])
    T = readtable('VSFP_50Hz_Proc.xlsx');
   
    % file saving destination
    FolderName = [startFile,'\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI]; 
    %check to see if Folder Exists to Save data, if not create it
    if ~exist(FolderName, 'dir')
      mkdir(FolderName)
    end
  
    all_mice = unique(T.Mouse); 
    all_mice(2) = []; % remove VSFP 25, no pupil data
   
    HV_All = [];
    PV_All = [];
    PH_All = [];
    FV_All = [];
    FH_All = [];
    muHVC1 = [];
    muPVC1 = [];
    muPHC1 = [];
    muFVC1 = [];
    muFHC1 = [];
    stdHVC1 = [];
    stdPVC1 = [];
    stdPHC1 = [];
    stdFVC1 = [];
    stdFHC1 = [];
    num = 0;
    
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        muHV = []; 
        muPV = []; 
        muPH = []; 
        muFV = []; 
        muFH = []; 
        
        %%% select an ROI from the Volt corr map for the data from
        %load in avg corr map generated from VSFP_CorrMap_Pupil.m, if it does not
        %exist will run that code first 
        disp('Checking to see if .mat Avg Corr Map already exists')
        dataFolder = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil\'];
        fname = [mouse,'_','avgCorrMap.mat'];
        if isfile(fullfile(dataFolder, fname))
           disp('File Exists!')
           load(fullfile(dataFolder, fname), 'avgTrial_Volt');
        else 
           disp('File does not exists, running VSFP_CorrMap_Pupil.m to generate needed file')
           VSFP_CorrMap_Pupil(mouse)
           load(fullfile(dataFolder, fname), 'avgTrial_Volt');
        end
        disp('Select desired ROI from cortical map')
        avgTrial_Volt = mean(avgTrial_Volt,3);
        imagesc(avgTrial_Volt), title('Select ROI: Avg Corr Map for Volt and Pupil'), colorbar
        r1 = drawcircle('Label','ROI','Color',[1 0 0]);
        center1 = r1. Center; 
        radius1 = r1.Radius;
        maskImg = createMask(r1);

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
            face_data2x = zscore(vsfp_data.face_data2x)'; 
            
            %voltage signal 
            data_volt = vsfp_data.out.imgDR3;

            %hemo signal 
            data_hemo = -1 * vsfp_data.out.hemoLP;

            % mask out all other areas
            % get correct pixel size
            maskImg = imresize(maskImg, [vsfp_data.out.sX, vsfp_data.out.sY]); 
            maskImg = maskImg(:,:,1);

            %apply mask on volt and hemo data 
            masked_volt = imrotate(data_volt, 270) .* maskImg;
            masked_hemo = imrotate(data_hemo, 270).* maskImg;

            % Find the regional "Global Signal"
            volt2D = reshape(masked_volt, 10000,[]);
            voltTemp = zscore(volt2D,[],2);
            gsVolt = mean(voltTemp,1);

            hemo2D = reshape(masked_hemo, 10000,[]);
            hemoTemp = zscore(hemo2D,[],2);
            gsHemo = mean(hemoTemp,1);

            Fs = vsfp_data.out.Fs;  % Sampling Frequency (Hz)
            t = 0:1/Fs:length(pupil_data2x)/Fs-1/Fs;

            % plot both signals in the ROI 
            f1 = figure(1);
            subplot(3,1,1)
            plot(t, gsHemo(1:length(pupil_data2x)),'DisplayName','Avg Global Hemo')
            hold on
            plot(t, gsVolt(1:length(pupil_data2x)), 'DisplayName','Avg Global Volt')
            legend()
            subplot(3,1,2)
            plot(t,pupil_data2x);
            xlabel('Time (s)')
            title('Pupil Diameter')
            subplot(3,1,3)
            plot(t(1:length(face_data2x)),face_data2x);
            xlabel('Time (s)')
            title('Orofacial Movements')
            
            f2 = figure(2);
            cwt(gsVolt(1:length(pupil_data2x)),Fs);
            colormap(color)
            title('Magnitude Scaliogram Volt')
            
            f3 = figure(3);
            cwt(gsHemo(1:length(pupil_data2x)),Fs);
            colormap(color)
            title('Magnitude Scaliogram Hemo')

            f4 = figure(4);
            cwt(pupil_data2x ,Fs);
            colormap(color)
            title('Pupil Diameter Wavelet')

            f5 = figure(5);
            cwt(face_data2x,Fs);
            colormap(color)
            title('Orofacial Movement Wavelet')
            
            f6 = figure(6);
            % take data from 20 to 60 ms (110 samples)
            s = find(t == 20);
            s1 = find( round(t, 2) == 60);
            subplot(6,1,1)
            [wcoh,~,F,coi] = wcoherence(gsHemo(1:length(pupil_data2x)),gsVolt(1:length(pupil_data2x)),Fs);
            helperPlotCoherence(wcoh(:,s:s1),t(1,s:s1),F,coi(s:s1),'Time (s)','Frequency (Hz)')
            title('Wavelet Coherence Hemo & Volt')
            muHV = [muHV,mean(wcoh(1:110,s:s1),2)];
            colormap(color)
            subplot(6,1,2)
            [wcoh,~,F,coi] = wcoherence(pupil_data2x,gsVolt(1:length(pupil_data2x)),Fs);
            helperPlotCoherence(wcoh(:,s:s1),t(1,s:s1),F,coi(s:s1),'Time (s)','Frequency (Hz)')
            title('Wavelet Coherence Pupil & Volt')
            muPV = [muPV,mean(wcoh(1:110,s:s1),2)];
            colormap(color)
            subplot(6,1,3)
            [wcoh,~,F,coi] = wcoherence(pupil_data2x,gsHemo(1:length(pupil_data2x)),Fs);
            helperPlotCoherence(wcoh(:,s:s1),t(1,s:s1),F,coi(s:s1),'Time (s)','Frequency (Hz)')
            title('Wavelet Coherence Pupil & Hemo')
            muPH = [muPH,mean(wcoh(1:110,s:s1),2)];
            colormap(color)
            
            % plot coherence for face data
            subplot(6,1,4)
            [wcoh,~,F,coi] = wcoherence(face_data2x,gsVolt(1:length(face_data2x)),Fs);
            helperPlotCoherence(wcoh(:,s:s1),t(1,s:s1),F,coi(s:s1),'Time (s)','Frequency (Hz)')
            title('Wavelet Coherence Face & Volt')
            muFV = [muFV,mean(wcoh(1:110,s:s1),2)];
            colormap(color)
            subplot(6,1,5)
            [wcoh,~,F,coi] = wcoherence(face_data2x,gsHemo(1:length(face_data2x)),Fs);
            helperPlotCoherence(wcoh(:,s:s1),t(1,s:s1),F,coi(s:s1),'Time (s)','Frequency (Hz)')
            title('Wavelet Coherence Face & Hemo')
            muFH = [muFH,mean(wcoh(1:110,s:s1),2)];
            colormap(color)
            % plot raw traces
            xlim([20 60])
            subplot(6,1,6)
            plot(pupil_data2x(s:s1),'Color', 'k','DisplayName','Pupil')
            hold on 
            plot(face_data2x(s:s1),'Color', 'k','DisplayName','Face')
            plot(gsVolt(s:s1),'Color', '#0072BD', 'DisplayName','Volt')
            plot(gsHemo(s:s1),'Color', '#A2142F', 'DisplayName','Hemo')
            title('Signal Traces')
            colorbar
            xlim([0 2000])
            legend
            hold off
            %axis equal
            
            %% save the made figures
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_ROI'];
            saveFig(f1, fname, FolderName);
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_Volt_Scalogram'];
            saveFig(f2, fname, FolderName);
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_Hemo_Scalogram'];
            saveFig(f3, fname, FolderName);
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_Pupil_Scalogram'];
            saveFig(f4, fname, FolderName);
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_Orofacial_Scalogram'];
            saveFig(f5, fname, FolderName);
            fname = [mouse,'_', num2str(date) ,'_',num2str(trial),'_Wavelet_Coherence'];
            saveFig(f6, fname, FolderName);

            close(f1,f2,f3,f4,f5,f6)
        end

        display(m);
        muHV_all = mean(muHV(32:110,:), 2);
        stdHV_all = (std(muHV(32:110,:)')./sqrt(size(muHV,2)));
        muPV_all = mean(muPV(32:110,:),2);
        stdPV_all = (std(muPV(32:110,:)')./sqrt(size(muPV,2)));
        muPH_all = mean(muPH(32:110,:),2);
        stdPH_all = (std(muPH(32:110,:)')./sqrt(size(muPH,2)));

        muFV_all = mean(muFV(32:110,:),2);
        stdFV_all = (std(muFV(32:110,:)')./sqrt(size(muFV,2)));
        muFH_all = mean(muFH(32:110,:),2);
        stdFH_all = (std(muFH(32:110,:)')./sqrt(size(muFH,2)));
        
        f7 = figure(7);
        subplot(5,1,1)
        for i = 1:size(muHV,2)
            plot(F(32:110), muHV(32:110,i), 'Color', [0.5 0.5 0.5])
            hold on
        end
        errorbar(F(32:110),mean(muHV(32:110,:),2), (std(muHV(32:110,:)')./sqrt(size(muHV,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
        title('Wavlet Coherence Hemo and Volt')
        xlabel('Freq')
        ylabel('Magnitude Squared Coherence') 
        
        subplot(5,1,2)
        for i = 1:size(muPV,2)
            plot(F(32:110), muPV(32:110,i), 'Color', [0.5 0.5 0.5])
            hold on
        end
        errorbar(F(32:110),mean(muPV(32:110,:),2), (std(muPV(32:110,:)')./sqrt(size(muPV,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
        title('Wavlet Coherence Pupil and Volt')
        xlabel('Freq')
        ylabel('Magnitude Squared Coherence') 
        
        subplot(5,1,3)
        for i = 1:size(muPH,2)
            plot(F(32:110), muPH(32:110,i), 'Color', [0.5 0.5 0.5])
            hold on
        end
        errorbar(F(32:110),mean(muPH(32:110,:),2), (std(muPH(32:110,:)')./sqrt(size(muPH,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
        title('Wavlet Coherence Pupil and Hemo')
        xlabel('Freq')
        ylabel('Magnitude Squared Coherence') 

        subplot(5,1,4)
        for i = 1:size(muFV,2)
            plot(F(32:110), muFV(32:110,i), 'Color', [0.5 0.5 0.5])
            hold on
        end
        errorbar(F(32:110),mean(muFV(32:110,:),2), (std(muFV(32:110,:)')./sqrt(size(muFV,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
        title('Wavlet Coherence Face and Volt')
        xlabel('Freq')
        ylabel('Magnitude Squared Coherence') 
        
        subplot(5,1,5)
        for i = 1:size(muFH,2)
            plot(F(32:110), muFH(32:110,i), 'Color', [0.5 0.5 0.5])
            hold on
        end
        errorbar(F(32:110),mean(muFH(32:110,:),2), (std(muFH(32:110,:)')./sqrt(size(muFH,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
        title('Wavlet Coherence Face and Hemo')
        xlabel('Freq')
        ylabel('Magnitude Squared Coherence') 
        hold off
       
        fname = ['Avg_', mouse,'_0', num2str(date),'_Wavelet_Coherence'];
        saveFig(f7, fname, FolderName);  
        close all; 

        %save variables to avg across all trials 
        HV_All = [HV_All, muHV];
        PV_All = [PV_All, muPV];
        PH_All = [PH_All, muPH];
        FV_All = [FV_All, muFV];
        FH_All = [FH_All, muFH];
        
        %save mu and std to plot at end to compare 
        muHVC1 = cat(2,muHVC1,muHV_all);
        stdHVC1 = cat(1,stdHVC1, stdHV_all);
        muPVC1 = cat(2,muPVC1,muPV_all);
        stdPVC1 = cat(1,stdPVC1, stdPV_all);
        muPHC1 = cat(2,muPHC1,muPH_all);
        stdPHC1 = cat(1,stdPHC1, stdPH_all);
        muFVC1 = cat(2,muFVC1,muFV_all);
        stdFVC1 = cat(1,stdFVC1, stdFV_all);
        muFHC1 = cat(2,muFHC1,muFH_all);
        stdFHC1 = cat(1,stdFHC1, stdFH_all);
        num = num + k;
    end

    close all;
    %save the files that contain data for each trial needed to significance
    %testing
    saveData = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI, 'Coherence_HV'];  % Your destination folder
    save(saveData, 'HV_All')
    saveData = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI,  'Coherence_PV'];  % Your destination folder
    save(saveData, 'PV_All')
    saveData = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI,  'Coherence_PH'];  % Your destination folder
    save(saveData, 'PH_All')
    saveData = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI,  'Coherence_FV'];  % Your destination folder
    save(saveData, 'FV_All')
    saveData = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\', ROI,  'Coherence_FH'];  % Your destination folder
    save(saveData, 'FH_All')
    
    f8 = figure(8);
    subplot(5,1,1)
    for i = 1:size(HV_All,2)
        plot(F(32:110), HV_All(32:110,i), 'Color', [0.7 0.7 0.7])
        hold on
    end
    err = std(HV_All(32:110,:)')./sqrt(size(HV_All,2));
    shade(F(32:110), mean(HV_All(32:110,:),2)+err',F(32:110),mean(HV_All(32:110,:),2)- err', 'FillType',[1 2])
 	hold on 
 	plot(F(32:110), mean(HV_All(32:110,:),2), 'k', 'LineWidth',2)
%    errorbar(F(32:110),mean(HV_All(32:110,:),2), (std(HV_All(32:110,:)')./sqrt(size(HV_All,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    title('Wavlet Coherence Hemo and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,2)
    for i = 1:size(PV_All,2)
        plot(F(32:110), PV_All(32:110,i), 'Color', [0.7 0.7 0.7])
        hold on
    end
    err = std(PV_All(32:110,:)')./sqrt(size(PV_All,2));
    shade(F(32:110), mean(PV_All(32:110,:),2)+err',F(32:110),mean(PV_All(32:110,:),2)- err', 'FillType',[1 2])
	hold on 
	plot(F(32:110), mean(PV_All(32:110,:),2), 'k', 'LineWidth',2)
%     errorbar(F(32:110),mean(PV_All(32:110,:),2), (std(PV_All(32:110,:)')./sqrt(size(PV_All,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    title('Wavlet Coherence Pupil and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,3)
    for i = 1:size(PH_All,2)
        plot(F(32:110), PH_All(32:110,i), 'Color', [0.7 0.7 0.7])
        hold on
    end
    err = std(PH_All(32:110,:)')./sqrt(size(PH_All,2));
    shade(F(32:110), mean(PH_All(32:110,:),2)+err',F(32:110),mean(PH_All(32:110,:),2)- err', 'FillType',[1 2])
	hold on 
	plot(F(32:110), mean(PH_All(32:110,:),2), 'k', 'LineWidth',2)
%     errorbar(F(32:110),mean(PH_All(32:110,:),2), (std(PH_All(32:110,:)')./sqrt(size(PH_All,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    title('Wavlet Coherence Pupil and Hemo')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,4)
    for i = 1:size(FH_All,2)
        plot(F(32:110), FH_All(32:110,i), 'Color', [0.7 0.7 0.7])
        hold on
    end
    err = std(FH_All(32:110,:)')./sqrt(size(FH_All,2));
    shade(F(32:110), mean(FH_All(32:110,:),2)+err',F(32:110),mean(FH_All(32:110,:),2)- err', 'FillType',[1 2])
	hold on 
	plot(F(32:110), mean(FH_All(32:110,:),2), 'k', 'LineWidth',2)
%     errorbar(F(32:110),mean(PH_All(32:110,:),2), (std(PH_All(32:110,:)')./sqrt(size(PH_All,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    title('Wavlet Coherence Face and Hemo')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,5)
    for i = 1:size(FV_All,2)
        plot(F(32:110), FV_All(32:110,i), 'Color', [0.7 0.7 0.7])
        hold on
    end
    err = std(FV_All(32:110,:)')./sqrt(size(FV_All,2));
    shade(F(32:110), mean(FV_All(32:110,:),2)+err',F(32:110),mean(FV_All(32:110,:),2)- err', 'FillType',[1 2])
	hold on 
	plot(F(32:110), mean(FV_All(32:110,:),2), 'k', 'LineWidth',2)
%     errorbar(F(32:110),mean(PH_All(32:110,:),2), (std(PH_All(32:110,:)')./sqrt(size(PH_All,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    title('Wavlet Coherence Face and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off
    
    

    % plot all the trial avgs with error bars
    f9 = figure(9);
    subplot(5,1,1)
    for i = 1:size(muHVC1,2)
        errorbar(F(32:110),muHVC1(:,i), stdHVC1(i,:), 'Color', '#7E2F8E', 'LineWidth', 1)
        hold on
    end

    title('Wavlet Coherence Hemo and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,2)
    for i = 1:size(muPVC1,2)
        errorbar(F(32:110),muPVC1(:,i), stdPVC1(i,:), 'Color', '#7E2F8E', 'LineWidth', 1)
        hold on 
    end
    title('Wavlet Coherence Pupil and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,3)
    for i = 1:size(muPHC1,2)
        errorbar(F(32:110), muPHC1(:,i), stdPHC1(i,:), 'Color', '#7E2F8E', 'LineWidth', 1)
        hold on 
    end
    title('Wavlet Coherence Pupil and Hemo')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    subplot(5,1,4)
    for i = 1:size(muFHC1,2)
        errorbar(F(32:110), muFHC1(:,i), stdFHC1(i,:), 'Color', '#7E2F8E', 'LineWidth', 1)
        hold on 
    end
    title('Wavlet Coherence Face and Hemo')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off
   
    subplot(5,1,5)
    for i = 1:size(muFVC1,2)
        errorbar(F(32:110), muFVC1(:,i), stdFVC1(i,:), 'Color', '#7E2F8E', 'LineWidth', 1)
        hold on 
    end
    title('Wavlet Coherence Face and Volt')
    xlabel('Freq')
    ylabel('Magnitude Squared Coherence') 
    hold off

    fname = ['All_','Wavelet_Coherence.fig'];
    saveFig(f9, fname,FolderName);
    fname = ['All_per_trial_','Wavelet_Coherence.fig'];
    saveFig(f8, fname,FolderName);
    
end


