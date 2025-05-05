% Function takes xlsx sheet and generates figures for following analysis
%   1) Runs significance test on imaging data
%       -- Can be done for Volt and Hemo versus shuffled control
%       -- Can be done between volt and hemo for LP and BP rest_mvmnt data
%       -- Looks at rest and mvmt versus shuffled control data 
%
% Inputs: 
%        shuffled -> runs test on shuffled controls versus true corr maps
%               VSFP_CorrMap_MannWhitneyTest('shuffled')
%
%               VSFP_CorrMap_MannWhitneyTest('rest_mvmt','all','volt')
%               VSFP_CorrMap_MannWhitneyTest('rest_mvmt','LP','hemo')
%               VSFP_CorrMap_MannWhitneyTest('rest_mvmt','BP','volt')
%
%               VSFP_CorrMap_MannWhitneyTest('rest_mvmt_shuffled','LP','hemo',)
%               VSFP_CorrMap_MannWhitneyTest('rest_mvmt_shuffled','BP','volt')
% Outputs: 
%        .fig files of analysis 
%        .pdf file of analysis
%        .mat files that contains corr values for 100 runs
%
% Written by Lisa Meyer-Baese

function VSFP_CorrMap_MannWhitneyTest(varargin)
    
    close all; 
    
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    load('YIOrBr8.mat')
    load('CustomColorMapPaper.mat')

    % select image used for mask
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    
    % load in cortical map Cortical_Map.png
    disp('Selecting Cortical Map to Use for Registration')
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));


    %% for shuffled 
    if isequal(varargin{1},'shuffled')
        
        % load shuffled data 
        addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil_shuffled'])
        %addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Face_shuffled'])
        load('AvgTrial_Volt.mat')
        shuffleVolt = avgTrial_Volt;
        load('AvgTrial_Hemo.mat')
        shuffleHemo = avgTrial_Hemo;
        
        %load real corr maps 
        allVolt = [];
        allHemo = [];
        
        addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil'])
        %addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Face'])
        load('VSFP24_avgCorrMap.mat', 'avgTrial_Volt')
        allVolt = cat(3,allVolt,avgTrial_Volt);
        load('VSFP24_avgCorrMap.mat', 'avgTrial_Hemo')
        allHemo = cat(3,allHemo,avgTrial_Hemo);
        
        load('VSFP27_avgCorrMap.mat', 'avgTrial_Volt')
        allVolt = cat(3,allVolt,avgTrial_Volt);
        load('VSFP27_avgCorrMap.mat', 'avgTrial_Hemo')
        allHemo = cat(3,allHemo,avgTrial_Hemo);
        
        load('VSFP28_avgCorrMap.mat', 'avgTrial_Volt')
        allVolt = cat(3,allVolt,avgTrial_Volt);
        load('VSFP28_avgCorrMap.mat', 'avgTrial_Hemo')
        allHemo = cat(3,allHemo,avgTrial_Hemo);
        
        load('VSFP29_avgCorrMap.mat', 'avgTrial_Volt')
        allVolt = cat(3,allVolt,avgTrial_Volt);
        load('VSFP29_avgCorrMap.mat', 'avgTrial_Hemo')
        allHemo = cat(3,allHemo,avgTrial_Hemo);
        
        load('VSFP30_avgCorrMap.mat', 'avgTrial_Volt')
        allVolt = cat(3,allVolt,avgTrial_Volt);
        load('VSFP30_avgCorrMap.mat', 'avgTrial_Hemo')
        allHemo = cat(3,allHemo,avgTrial_Hemo);
    
        % make map for voltage data 
        pValMap = zeros(100,100);
        hValMap = zeros(100,100);
        bfCorrected = 0.05/10000;
        for x = 1:100 
            for y = 1:100 
                shuffle = squeeze(shuffleVolt(x,y,:));
                real = squeeze(allVolt(x,y,:));
                [p,h] = ranksum(shuffle,real, 'alpha',bfCorrected, 'tail','both');
                pValMap(x,y) = p;
                hValMap(x,y) = h;
            end
        end
    
        % mask out non-cortex
        load('YIOrBr8.mat')
        % get correct pixel size
        maskImg = imresize(maskImg, [100,100]); 
        maskImg = maskImg(:,:,1);

        %apply mask on transformed data 
        maskk = maskImg < 100;
        pValMap(maskk)= 0;
        hValMap(maskk) = 0;

        % get atlas to overlay
        fixed_map = imresize(fixed_map, [100, 100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
    
        bottom = -1;
        top = 1;
        f1 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,3,1)
        imagesc(mean(allVolt,3)); colorbar, caxis([-0.1,0.6]), colormap(CustomColormap1)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Volt')
        axis square
        axis off
    
        subplot(1,3,2)
        imagesc(mean(shuffleVolt,3)), colorbar, caxis([-0.1,0.6]), colormap(CustomColormap1)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of Shuffled Control')
        axis square
        axis off

        %find the difference between volt and shuffled to split color scale on p
        %value map based on + and - differences 
        diff = double(mean(allVolt,3) > mean(shuffleVolt,3));
        diff(diff == 0) = -1;
        newP = 0.05 / 100;
        ind = double(pValMap <= newP);
        subplot(1,3,3)
        imagesc(ind.*diff), colorbar, caxis([-1,1])
        colormap(CustomColormap1)    
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values')
        axis square
        axis off
    

        % make map for hemo data 
        pValMap = zeros(100,100);
        hValMap = zeros(100,100);
        for x = 1:100 
            for y = 1:100 
                shuffle = squeeze(shuffleHemo(x,y,:));
                real = squeeze(allHemo(x,y,:));
                [p,h] = ranksum(shuffle,real, 'alpha',0.05, 'tail','both');
                pValMap(x,y) = p;
                hValMap(x,y) = h;
            end
        end

        % mask out non-cortex
        % get correct pixel size
        maskImg = imresize(maskImg, [100,100]); 
        maskImg = maskImg(:,:,1);

        %apply mask on transformed data 
        maskk = maskImg < 100;
        pValMap(maskk)= 0;
        hValMap(maskk) = 0;

        % get atlas to overlay
        fixed_map = imresize(fixed_map, [100, 100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts

        bottom = -1;
        top = 1;
        f2 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(1,3,1)
        imagesc(mean(allHemo,3)), colorbar, caxis([-0.5,0.3]), colormap(CustomColormap1)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Hemo')
        axis square
        axis off

        subplot(1,3,2)
        imagesc(mean(shuffleHemo,3)),  colorbar, caxis([-0.5,0.3]),colormap(CustomColormap1)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of Shuffled Control')
        axis square
        axis off

        %find the difference between volt and shuffled to split color scale on p
        %value map based on + and - differences 
        diff = double(mean(allHemo,3) > mean(shuffleHemo,3));
        diff(diff == 0) = -1;
        newP = 0.05 / 100;
        ind = double(pValMap <= newP);
        subplot(1,3,3)
        imagesc(ind.*diff), colorbar, caxis([-1,1])
        colormap(CustomColormap1)    
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values')
        axis square
        axis off

        FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil_shuffled\'];  % Your destination folder
        %FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Face_shuffled\'];  % Your destination folder
        fname = 'Volt_CorrMap_StatTest';
        saveFig(f1, fname, FolderName);
        fname = 'Hemo_CorrMap_StatTest';
        saveFig(f2, fname, FolderName);
        
    %% for mvmt vs rest    
    elseif isequal(varargin{1},'rest_mvmt')
        mvmt_all = [];
        rest_all = [];
        if isequal(varargin{2},'all')
            addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Face\2 sec\Avg Maps'])
            if isequal(varargin{3},'volt')
                load('VSFP24_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP24_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP27_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP27_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP28_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP28_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP29_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP29_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP30_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP30_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);
            elseif isequal(varargin{3},'hemo')       
                load('VSFP24_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP24_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP27_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP27_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP28_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP28_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP29_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP29_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP30_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP30_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);
            end
        % If Low Pass Filter Data 
        elseif isequal(varargin{2},'LP')
            arr = [-0.2, 0.6];
            addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\LP_2\'])
            if isequal(varargin{3},'volt')
                load('VSFP24_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP24_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP27_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP27_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP28_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP28_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP29_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP29_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP30_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP30_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);
            elseif isequal(varargin{3},'hemo')       
                load('VSFP24_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP24_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP27_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP27_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP28_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP28_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP29_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP29_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP30_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP30_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);
            end
        elseif isequal(varargin{2},'BP')
            arr = [-0.1, 0.3];
             addpath([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\BP_2\'])
            if isequal(varargin{3},'volt')
                load('VSFP24_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP24_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP27_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP27_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP28_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP28_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP29_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP29_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP30_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP30_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);
            elseif isequal(varargin{3},'hemo')       
                load('VSFP24_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP24_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP27_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP27_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP28_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP28_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP29_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP29_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP30_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP30_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);
            end
        end
        
        % make map for voltage data 
        pValMap = zeros(100,100);
        hValMap = zeros(100,100);
        for x = 1:100 
            for y = 1:100 
                shuffle = squeeze(mvmt_all(x,y,:));
                real = squeeze(rest_all(x,y,:));
                [p,h] = ranksum(shuffle,real, 'alpha',0.05, 'tail','both');
                pValMap(x,y) = p;
                hValMap(x,y) = h;
            end
        end
    
       % mask out non-cortex
        % get correct pixel size
        maskImg = imresize(maskImg, [100,100]); 
        maskImg = maskImg(:,:,1);

        %apply mask on transformed data 
        maskk = maskImg < 100;
        pValMap(maskk)= 0;
        hValMap(maskk) = 0;

        % get atlas to overlay
        fixed_map = imresize(fixed_map, [100, 100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
    
        bottom = -1;
        top = 1;
        f1 = figure(1);
        subplot(2,2,1)
        imagesc(mean(rest_all,3)); colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Rest')
        axis square
        axis off
    
        subplot(2,2,2)
        imagesc(mean(mvmt_all,3)), colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Mvmt')
        axis square
        axis off

        subplot(2,2,3)
        imagesc(pValMap),  colorbar
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values')
        axis square
        axis off
    
        subplot(2,2,4)
        imagesc(hValMap),  colorbar
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('H Values')
        axis square
        axis off
         %find the difference between mvmt and rest to split color scale on p
        %value map based on + and - differences 
        diff = double(mean(mvmt_all,3) > mean(rest_all,3));
        diff(diff == 0) = -1;
        imagesc(diff)
        newP = 0.05 / 100;
        f2 = figure(2);
        ind = double(pValMap <= newP);
        imagesc(ind.*diff), colorbar, caxis([-1,1])
        colormap(CustomColormap1)    
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values')
        axis square
        axis off

        FolderName = ([startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\']);
        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3}];
        saveFig(f1, fname, FolderName);
        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3},'_p-value'];
        saveFig(f2, fname, FolderName);
        close all
    %% for mvmt vs rest  vs shuffled 
    elseif isequal(varargin{1},'rest_mvmt_shuffled')
        mvmt_all = [];
        rest_all = [];

        % If Low Pass Filter Data 
        if isequal(varargin{2},'LP')
            arr = [-0.2, 0.6];
            addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\LP_2'])
            if isequal(varargin{3},'volt')
                load('VSFP24_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP24_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP27_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP27_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP28_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP28_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP29_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP29_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP30_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP30_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);
                
                addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq_shuffled\LP_2'])
                load('avg_CorrMap_Volt_mvmt.mat')
                mvmt_shuffle = avgTrial_Volt_mvmt;
                load('avg_CorrMap_Volt_rest.mat')
                rest_shuffle = avgTrial_Volt_mvmt;
                
            elseif isequal(varargin{3},'hemo')       
                load('VSFP24_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP24_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP27_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP27_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP28_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP28_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP29_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP29_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP30_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP30_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);
                
                addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq_shuffled\LP_2'])
                load('avg_CorrMap_Hemo_mvmt.mat')
                mvmt_shuffle = avgTrial_Hemo_mvmt;
                load('avg_CorrMap_Hemo_rest.mat')
                rest_shuffle = avgTrial_Hemo_mvmt;
            end
        elseif isequal(varargin{2},'BP')
            arr = [-0.1, 0.3];
            addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq\BP_2'])
            if isequal(varargin{3},'volt')
                load('VSFP24_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP24_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP27_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP27_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP28_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP28_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP29_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP29_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);

                load('VSFP30_avg_CorrMap_Volt_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Volt_mvmt);
                load('VSFP30_avg_CorrMap_Volt_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Volt_rest);
                
                addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq_shuffled\BP_2'])
                load('avg_CorrMap_Volt_mvmt.mat')
                mvmt_shuffle = avgTrial_Volt_mvmt;
                load('avg_CorrMap_Volt_rest.mat')
                rest_shuffle = avgTrial_Volt_mvmt;
                
            elseif isequal(varargin{3},'hemo')       
                load('VSFP24_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP24_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP27_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP27_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP28_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP28_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP29_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP29_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);

                load('VSFP30_avg_CorrMap_Hemo_mvmt.mat')
                mvmt_all = cat(3, mvmt_all, avgTrial_Hemo_mvmt);
                load('VSFP30_avg_CorrMap_Hemo_rest.mat')
                rest_all = cat(3, rest_all, avgTrial_Hemo_rest);
                
                addpath ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq_shuffled\BP_2'])
                load('avg_CorrMap_Hemo_mvmt.mat')
                mvmt_shuffle = avgTrial_Hemo_mvmt;
                load('avg_CorrMap_Hemo_rest.mat')
                rest_shuffle = avgTrial_Hemo_mvmt;
            end
        end
        
        % make map for mvmt data 
        mvmt_pValMap = zeros(100,100);
        mvmt_hValMap = zeros(100,100);
        for x = 1:100 
            for y = 1:100 
                shuffle = squeeze(mvmt_shuffle(x,y,2:end));
                real = squeeze(mvmt_all(x,y,:));
                [p,h] = ranksum(shuffle,real, 'alpha',0.05, 'tail','both');
                mvmt_pValMap(x,y) = p;
                mvmt_hValMap(x,y) = h;
            end
        end
    
       % mask out non-cortex
        % get correct pixel size
        maskImg = imresize(maskImg, [100,100]); 
        maskImg = maskImg(:,:,1);

        %apply mask on transformed data 
        maskk = maskImg < 100;
        mvmt_pValMap(maskk)= 0;
        mvmt_hValMap(maskk) = 0;
        
        
        % make map for rest data 
        rest_pValMap = zeros(100,100);
        rest_hValMap = zeros(100,100);
        for x = 1:100 
            for y = 1:100 
                shuffle = squeeze(rest_shuffle(x,y,2:end));
                real = squeeze(rest_all(x,y,:));
                [p,h] = ranksum(shuffle,real, 'alpha',0.05, 'tail','both');
                rest_pValMap(x,y) = p;
                rest_hValMap(x,y) = h;
            end
        end
    
       % mask out non-cortex
        % get correct pixel size
        maskImg = imresize(maskImg, [100,100]); 
        maskImg = maskImg(:,:,1);

        %apply mask on transformed data 
        maskk = maskImg < 100;
        rest_pValMap(maskk)= 0;
        rest_hValMap(maskk) = 0;

        % get atlas to overlay
        fixed_map = imresize(fixed_map, [100, 100]);
        M = repmat(all(~fixed_map,3),[1 1 3]); %mask black parts
    
        
        bottom = -1;
        top = 1;
        f1 = figure(1);
        subplot(2,2,1)
        imagesc(mean(rest_shuffle,3)); colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Rest Shuffle')
        axis square
        axis off
    
        subplot(2,2,2)
        imagesc(mean(rest_all,3)), colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Rest')
        axis square
        axis off
        
        subplot(2,2,3)
        imagesc(mean(mvmt_shuffle,3)); colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Mvmt Shuffle')
        axis square
        axis off
    
        subplot(2,2,4)
        imagesc(mean(mvmt_all,3)), colorbar, caxis(arr)
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('Avg of All Mvmt')
        axis square
        axis off
        
        f2 = figure(2);
        
        subplot(2,1,1)
        %find the difference between mvmt and rest to split color scale on p
        %value map based on + and - differences 
        diff = double(mean(rest_all,3) > mean(rest_shuffle,3));
        diff(diff == 0) = -1;
        imagesc(diff)
        newP = 0.05 / 100;
        ind = double(rest_pValMap <= newP);
        imagesc(ind.*diff), colorbar, caxis([-1,1])
        colormap(CustomColormap1)    
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values Rest')
        axis square
        axis off
        
        subplot(2,1,2)
        %find the difference between mvmt and rest to split color scale on p
        %value map based on + and - differences 
        diff = double(mean(mvmt_all,3) > mean(mvmt_shuffle,3));
        diff(diff == 0) = -1;
        imagesc(diff)
        newP = 0.05 / 100;
        ind = double(mvmt_pValMap <= newP);
        imagesc(ind.*diff), colorbar, caxis([-1,1])
        colormap(CustomColormap1)    
        hold on
        a = imagesc(M);
        alpha(a,.2)
        title('P Values Mvmt')
        axis square
        axis off

        FolderName = ([startFile ,'\VSFP ButterFly\Data\VSFP_CorrMap_BehavioralState_Freq_shuffled']);  % Your destination folder
        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3},'.fig'];
        savefig(f1, fullfile(FolderName, fname));
        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3},'.png'];
        saveas(f1, fullfile(FolderName, fname));

        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3},'_p-value.fig'];
        savefig(f2, fullfile(FolderName, fname));
        fname = [varargin{1}, '_',varargin{2}, '_',varargin{3},'_p-value.png'];
        saveas(f2, fullfile(FolderName, fname));
        close all
   
    end
   
end 
    
