
function preProcEphys_Final
    
    %% Edits made to Miao's original script to load LFP and sync the data with imaging 
    % 1. Form the metalist or load metalist file
    % 2. for each set of matched data
        %a. load imaging data
        %b. load LFP data and cropped it based on DIO6 input
    % 3. sort the LFP data based on the physical depth of channels 
    % 4. Saves the sorted LFP data into the same file as the pre processed
    %    imaging file 
    % Do this in a different script --> Calculates the correlation between imaging data and LFP data
    
   

    %% Set up parameters
    
    LFP_fs = 20000; %acquisition rate of LFP,
    steepness = 0.95; % for filtering data
    duration = 20.48; %imaging trial duration
    LFP_sampleRate = 20000; % /s, for calculating x axis for plotting
    IMG_sampleRate = 200; % /s, for calculating x axis for plotting
    
    
    % use the function to match intan channel map with acute probe numbers
    % probe: 32 channel Cambridge Neurotech H8b
    % depthInd = 1 -> shallow; depthInd = 32 -> deeps
    
    [depthInd] = convert_IntanChannel_acuteProbe;
     %% Form the metalist or load metalist file
    % Or use existing metalis
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end
    
    mouseID = {'VSFPE3'};
    %mouseID = {'VSFPE3', 'VSFPE4'};

    addpath ([startFile, '\VSFP ButterFly\VSFP_LFP\VSFP_LFP_Data_Scripts'])
    addpath ([startFile, '\VSFP ButterFly\Code\PreProc'])
    %load('VSFPE3_1216_LFP_img_metalist.mat')
    %load('VSFPE4_1216_LFP_img_metalist.mat','metaList');
    
    for m = 1:length(mouseID)
        load([mouseID{m},'_1216_LFP_img_metalist.mat']);

        % Read files
        for trialN = 1:length(metaList)
        %for trialN = 37:46
            display(trialN)
    
            nameSave = metaList(trialN).imgA_name(1:6);
            if isequal(nameSave, 'VSFPE3')
                FolderPath = [startFile, '\VSFP ButterFly\VSFP_LFP\LFP\VSFPE3_1216_2021'];
            else
                FolderPath = [startFile,'\VSFP ButterFly\VSFP_LFP\LFP\VSFPE4_1216_2021'];
            end
            
            % read_Intan_RHD2000_file
            [board_dig_in_data, amplifier_data] = read_Intan_RHD2000_file(metaList(trialN).LFP_name, FolderPath);
            [depthInd] = convert_IntanChannel_acuteProbe;
            % get the mouse name, trial session and date
            temp = split(metaList(trialN).imgA_name,'_');
            fpre = temp{1};
            temp2 = (split(temp{2}, '-'));
            fDate = temp2{1};
            fNum = temp2{2};
    
            % onset and offset of the LFP data based on the digital input & imaging duration        
            %load DIO6 trace
            dig_in_trialTrigger = board_dig_in_data(2, :); 
            %find the onset index of the first TTL pulse
            onsetInd = find(dig_in_trialTrigger == 1,  1, 'first'); 
            %offset index is determined by the duration of the trial
            offsetInd = onsetInd + duration * LFP_sampleRate - 1; 
            LFP_data = amplifier_data(:, onsetInd:offsetInd)';
    
            % Sort LFP data by depth
            sortedLFP = zeros(size(LFP_data));
            for ind = 1:length(depthInd)
                lfp = LFP_data(:, ind);
                % sort by depth from shallow to deep (e.g. 1 is shallow and 32 is
                % deep)
                sortedLFP(:, depthInd(ind)) = lfp;
            end
    
            %Save LFP Data to new field in structure 
            metaList(trialN).LFPData = sortedLFP;
    
            %% Load and process the VSFP Imaging Data too if it hasn't been done so already
            out = preProcVSFP10L_Final(fpre, fDate, fNum, fpre);
            saveFilePath = [startFile,'\VSFP ButterFly\Data\VSFP_Ephys\'];
            fileName = [fpre,'_',fNum,'_', fDate];
            save([saveFilePath,fileName,'.mat'],'out', 'sortedLFP','-v7.3')
    
        end
    end
   
end

function [depthInd] = convert_IntanChannel_acuteProbe

% Script to get to the depthInd in order to sort the channels by depth
% Two conversions are needed before the channels are correctly sorted
%     1. electrode to headstage
%     2. headstage to intan
% from shallow to deep
orderDepth = [17, 16, 21, 12, 22, 8, 25, 6, 27, 2, 29, 3, 28, 5, 30, 9, ...
    32, 11, 31, 7, 26, 1, 24, 4, 20, 13, 19, 14, 23, 10 ,18, 15];

%adapter pin numbers
adapter_map = [12, 10, 8, 6, 2, 3, 7, 11, 32, 31, 29, 27, 25, 22, 23, 21, ...
    16, 15, 14, 13, 4, 1, 5, 9, 30, 28, 26, 24, 20, 19, 18, 17];
%headstage pin numbers, representing the recorded channel number of Intan
headstage_map = [23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, ...
    24, 25, 26, 27, 28, 29, 30, 31, 0, 1, 2, 3, 4, 5, 6, 7];

% e.g. in23 (from Inatan) is electrode #12

% change the headstage from [0, 31] to [1, 32]
headstage_map = headstage_map + 1;
converted_map = zeros(1, length(headstage_map));
depthInd = zeros(1, length(headstage_map));

% from [1:32], find the index of in## with the headstage map
% 
for ind = 1:length(headstage_map)
    % from [1:32], find the index of in## with the headstage map
    headstage_ind = find(headstage_map == ind);
    % because the physical location of the adapter pin and headstage pin
    % matches
    % the same headstage_ind works for adapter_map
    converted_map(1, ind) = adapter_map(headstage_ind);
    % find depth of in## (shallow to deep)
    % using this depthInd to sort the channels results in
    % 1 -> shallow; 32 -> deep;
    depthInd(1, ind) = find(orderDepth == converted_map(1, ind));
end
end