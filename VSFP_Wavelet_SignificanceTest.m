% Function takes .mat files generated from 
%   VSFP_Wavelet_Shuffle_100.m and VSFP_Wavelet.m
%   1) takes data and plots mean coherence values for both 
%   2) determines the distribution of the data and does correct
%   significance test for it 
%
% Inputs: 
%        n/a
%
% Outputs: 
%        .fig files of analysis 
%        .png file of analysis
%
% Written by Lisa Meyer-Baese

function VSFP_Wavelet_SignificanceTest(varargin)
    
    % load in data for 100 shuffled controls 
    close all;
    
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    %% load in shuffled contols 
    Folder = [startFile, ' \VSFP ButterFly\Data\VSFP_Wavelet_Shuffle'];
    cd (Folder)

    % for ROI 1 
    load('1.mat')
    muPH_Shuffled1 = muPH_1;
    muPV_Shuffled1 = muPV_1;
    muVH_Shuffled1 = muHV_1;
    muFV_Shuffled1 = muFV_1;
    muFH_Shuffled1 = muFV_1;
    

    %shuffled control ROI 2
    load('2.mat')
    muPH_Shuffled2 = muPH_2;
    muPV_Shuffled2 = muPV_2;
    muVH_Shuffled2 = muHV_2;
    muFV_Shuffled2 = muFV_2;
    muFH_Shuffled2 = muFV_2;
    
    %shuffled control ROI 3
    load('3.mat')
    muPH_Shuffled3 = muPH_3;
    muPV_Shuffled3 = muPV_3;
    muVH_Shuffled3 = muHV_3;
    muFV_Shuffled3 = muFV_3;
    muFH_Shuffled3 = muFV_3;
    
    %shuffled control ROI 4
    load('4.mat')
    muPH_Shuffled4 = muPH_4;
    muPV_Shuffled4 = muPV_4;
    muVH_Shuffled4 = muHV_4;
    muFV_Shuffled4 = muFV_4;
    muFH_Shuffled4 = muFV_4;
    
    %% load in the VSFP coherence data for all trials 
    Folder = [startFile, ' \VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI'];
    cd(Folder)
    load('F.mat')
    
    % load in VSFP Coherence data for ROI 1
    Folder = [startFile,' \VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\VSFP_WaveletROI_one'];
    cd(Folder)
    load('Coherence_PH.mat')
    muPH_ROI1 = PH_All;
    load('Coherence_PV.mat')
    muPV_ROI1 = PV_All;
    load('Coherence_FH.mat')
    muFH_ROI1 = FH_All;
    load('Coherence_FV.mat')
    muFV_ROI1 = FV_All;
    load('Coherence_HV.mat')
    muHV_ROI1 = HV_All;
    %cROI1 = '#EDB120';
    cROI1 = "#A2142F";
    
    % load in VSFP Coherence data for ROI 2
    Folder = [startFile, '\VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\VSFP_WaveletROI_two'];
    cd(Folder)
    load('Coherence_PH.mat')
    muPH_ROI2 = PH_All;
    load('Coherence_PV.mat')
    muPV_ROI2 = PV_All;
     load('Coherence_FH.mat')
    muFH_ROI2 = FH_All;
    load('Coherence_FV.mat')
    muFV_ROI2 = FV_All;
    load('Coherence_HV.mat')
    muHV_ROI2 = HV_All;
    %cROI2 = '#77AC30';
    cROI2 = "#EDB120";
    
    % load in VSFP Coherence data for ROI 3
    Folder = [startFile, ' \VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\VSFP_WaveletROI_three'];
    cd(Folder)
    load('Coherence_PH.mat')
    muPH_ROI3 = PH_All;
    load('Coherence_PV.mat')
    muPV_ROI3 = PV_All;
    load('Coherence_FH.mat')
    muFH_ROI3 = FH_All;
    load('Coherence_FV.mat')
    muFV_ROI3 = FV_All;
    load('Coherence_HV.mat')
    muHV_ROI3 = HV_All;
    %cROI3 = '#0072BD';
    cROI3 = "#0072BD";
    
    % load in VSFP Coherence data for ROI 4
    Folder = [startFile,' \VSFP ButterFly\Data\VSFP_WaveletCoherence_ROI\VSFP_WaveletROI_four'];
    cd(Folder)
    load('Coherence_PH.mat')
    muPH_ROI4 = PH_All;
    load('Coherence_PV.mat')
    muPV_ROI4 = PV_All;
    load('Coherence_FH.mat')
    muFH_ROI4 = FH_All;
    load('Coherence_FV.mat')
    muFV_ROI4 = FV_All;
    load('Coherence_HV.mat')
    muHV_ROI4 = HV_All;
    %cROI4 = '#000000';
    cROI4 = "#77AC30";
    
    cd ([startFile,'\VSFP ButterFly\Code\Analysis'])
    %% plot the mean coherence for the two data sets
    %range8Hz = (1:110); %all freq
    %range8Hz = (51:110);  %up to 1Hz --> to 51
    %range8Hz = (11:24); % getting only 5-10Hz
    range8Hz = (13:110);  % up to 9 Hz --> (13:110) 
    range15Hz = (51:110);  %up to 1Hz --> to 51
    %range15Hz = (1:110); %all freq
    
    % plot the avg for all the 3 coherence pairs 
    f1 = figure(1);
    subplot(5,1,1)
    %errorbar(F(range8Hz),mean(muHV(range8Hz,:),2), (std(muHV(range8Hz,:)')./sqrt(size(muHV,2))), 'Color', '#7E2F8E', 'LineWidth', 1)
    hold on 
    %errorbar(F(range8Hz),mean(muVH_Shuffled(range8Hz,:),2), (std(muVH_Shuffled(range8Hz,:)')./sqrt(size(muVH_Shuffled,2))), 'Color', '#D95319', 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muHV_ROI1(range8Hz,:),2), (std(muHV_ROI1(range8Hz,:)')./sqrt(size(muHV_ROI1,2))), 'Color', cROI1, 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muHV_ROI2(range8Hz,:),2), (std(muHV_ROI2(range8Hz,:)')./sqrt(size(muHV_ROI2,2))), 'Color', cROI2, 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muHV_ROI3(range8Hz,:),2), (std(muHV_ROI3(range8Hz,:)')./sqrt(size(muHV_ROI3,2))), 'Color', cROI3, 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muHV_ROI4(range8Hz,:),2), (std(muHV_ROI4(range8Hz,:)')./sqrt(size(muHV_ROI4,2))), 'Color', cROI4, 'LineWidth', 1)
    
    errorbar(F(range8Hz),mean(muVH_Shuffled1(range8Hz,:),2), (std(muVH_Shuffled1(range8Hz,:)')./sqrt(size(muVH_Shuffled1,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muVH_Shuffled2(range8Hz,:),2), (std(muVH_Shuffled2(range8Hz,:)')./sqrt(size(muVH_Shuffled2,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muVH_Shuffled3(range8Hz,:),2), (std(muVH_Shuffled3(range8Hz,:)')./sqrt(size(muVH_Shuffled3,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range8Hz),mean(muVH_Shuffled4(range8Hz,:),2), (std(muVH_Shuffled4(range8Hz,:)')./sqrt(size(muVH_Shuffled4,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    title('Wavlet Coherence Hemo and Volt')
    xlabel('Freq')
    ylabel('Mag. Sq. Coherence')  

    subplot(5,1,2)
    hold on
    errorbar(F(range15Hz),mean(muPV_ROI1(range15Hz,:),2), (std(muPV_ROI1(range15Hz,:)')./sqrt(size(muPV_ROI1,2))), 'Color', cROI1, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_ROI2(range15Hz,:),2), (std(muPV_ROI2(range15Hz,:)')./sqrt(size(muPV_ROI2,2))), 'Color', cROI2, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_ROI3(range15Hz,:),2), (std(muPV_ROI3(range15Hz,:)')./sqrt(size(muPV_ROI3,2))), 'Color', cROI3, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_ROI4(range15Hz,:),2), (std(muPV_ROI4(range15Hz,:)')./sqrt(size(muPV_ROI4,2))), 'Color', cROI4, 'LineWidth', 1)
    
    errorbar(F(range15Hz),mean(muPV_Shuffled1(range15Hz,:),2), (std(muPV_Shuffled1(range15Hz,:)')./sqrt(size(muPV_Shuffled1,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_Shuffled2(range15Hz,:),2), (std(muPV_Shuffled2(range15Hz,:)')./sqrt(size(muPV_Shuffled2,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_Shuffled3(range15Hz,:),2), (std(muPV_Shuffled3(range15Hz,:)')./sqrt(size(muPV_Shuffled3,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPV_Shuffled4(range15Hz,:),2), (std(muPV_Shuffled4(range15Hz,:)')./sqrt(size(muPV_Shuffled4,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    title('Wavlet Coherence Pupil and Volt')
    xlabel('Freq')
    ylabel('Mag. Sq. Coherence') 

    subplot(5,1,3)
    hold on
    errorbar(F(range15Hz),mean(muPH_ROI1(range15Hz,:),2), (std(muPH_ROI1(range15Hz,:)')./sqrt(size(muPH_ROI1,2))),'Color', cROI1, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_ROI2(range15Hz,:),2), (std(muPH_ROI2(range15Hz,:)')./sqrt(size(muPH_ROI2,2))), 'Color',cROI2,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_ROI3(range15Hz,:),2), (std(muPH_ROI3(range15Hz,:)')./sqrt(size(muPH_ROI3,2))),'Color', cROI3,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_ROI4(range15Hz,:),2), (std(muPH_ROI4(range15Hz,:)')./sqrt(size(muPH_ROI4,2))),'Color', cROI4, 'LineWidth', 1)
    
    errorbar(F(range15Hz),mean(muPH_Shuffled1(range15Hz,:),2), (std(muPH_Shuffled1(range15Hz,:)')./sqrt(size(muPH_Shuffled1,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_Shuffled2(range15Hz,:),2), (std(muPH_Shuffled2(range15Hz,:)')./sqrt(size(muPH_Shuffled2,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_Shuffled3(range15Hz,:),2), (std(muPH_Shuffled3(range15Hz,:)')./sqrt(size(muPH_Shuffled3,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muPH_Shuffled4(range15Hz,:),2), (std(muPH_Shuffled4(range15Hz,:)')./sqrt(size(muPH_Shuffled4,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    title('Wavlet Coherence Pupil and Hemo')
    xlabel('Freq')
    ylabel('Mag. Sq. Coherence') 

    subplot(5,1,4)
    hold on
    errorbar(F(range15Hz),mean(muFH_ROI1(range15Hz,:),2), (std(muFH_ROI1(range15Hz,:)')./sqrt(size(muFH_ROI1,2))),'Color', cROI1, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_ROI2(range15Hz,:),2), (std(muFH_ROI2(range15Hz,:)')./sqrt(size(muFH_ROI2,2))), 'Color',cROI2,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_ROI3(range15Hz,:),2), (std(muFH_ROI3(range15Hz,:)')./sqrt(size(muFH_ROI3,2))),'Color', cROI3,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_ROI4(range15Hz,:),2), (std(muFH_ROI4(range15Hz,:)')./sqrt(size(muFH_ROI4,2))),'Color', cROI4, 'LineWidth', 1)
    
    errorbar(F(range15Hz),mean(muFH_Shuffled1(range15Hz,:),2), (std(muFH_Shuffled1(range15Hz,:)')./sqrt(size(muFH_Shuffled1,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_Shuffled2(range15Hz,:),2), (std(muFH_Shuffled2(range15Hz,:)')./sqrt(size(muFH_Shuffled2,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_Shuffled3(range15Hz,:),2), (std(muFH_Shuffled3(range15Hz,:)')./sqrt(size(muFH_Shuffled3,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFH_Shuffled4(range15Hz,:),2), (std(muFH_Shuffled4(range15Hz,:)')./sqrt(size(muFH_Shuffled4,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    title('Wavlet Coherence Face and Hemo')
    xlabel('Freq')
    ylabel('Mag. Sq. Coherence') 

    subplot(5,1,5)
    hold on
    errorbar(F(range15Hz),mean(muFV_ROI1(range15Hz,:),2), (std(muFV_ROI1(range15Hz,:)')./sqrt(size(muFV_ROI1,2))),'Color', cROI1, 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_ROI2(range15Hz,:),2), (std(muFV_ROI2(range15Hz,:)')./sqrt(size(muFV_ROI2,2))), 'Color',cROI2,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_ROI3(range15Hz,:),2), (std(muFV_ROI3(range15Hz,:)')./sqrt(size(muFV_ROI3,2))),'Color', cROI3,  'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_ROI4(range15Hz,:),2), (std(muFV_ROI4(range15Hz,:)')./sqrt(size(muFV_ROI4,2))),'Color', cROI4, 'LineWidth', 1)
    
    errorbar(F(range15Hz),mean(muFV_Shuffled1(range15Hz,:),2), (std(muFV_Shuffled1(range15Hz,:)')./sqrt(size(muFV_Shuffled1,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_Shuffled2(range15Hz,:),2), (std(muFV_Shuffled2(range15Hz,:)')./sqrt(size(muFV_Shuffled2,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_Shuffled3(range15Hz,:),2), (std(muFV_Shuffled3(range15Hz,:)')./sqrt(size(muFV_Shuffled3,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    errorbar(F(range15Hz),mean(muFV_Shuffled4(range15Hz,:),2), (std(muFV_Shuffled4(range15Hz,:)')./sqrt(size(muFV_Shuffled4,2))), 'Color', [.5 .5 .5], 'LineWidth', 1)
    title('Wavlet Coherence Face and Volt')
    xlabel('Freq')
    ylabel('Mag. Sq. Coherence') 
    
    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_Wavelet_SignificanceTest\'];
    fname = 'Avg Wavelet_Coherence';
    saveFig(f1,fname, FolderName);
    close(f1)
end