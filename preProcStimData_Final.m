% this function preprocesses all the resting state and whisker stimulation
% data collected on 8/22/23 for VSFP00 and VSFP01

function preProcStimData_Final
    %mice = {'VSFP00', 'VSFP01'};
    mice = {'VSFP01'};
    %BPOD_Files = {'VSFP01_Imaging_Whisker_Stim_Sweep_Frequency_20230918_140332.mat','VSFP01_Imaging_Whisker_Stim_Sweep_Frequency_20230918_144009.mat'};
    %BPOD_Files = {'VSFP00_Imaging_Whisker_Stim_Sweep_Frequency_20230917_200415.mat', 'VSFP00_Imaging_Whisker_Stim_Sweep_Frequency_20230917_204048.mat'};
    trials = [1:100];
    %side = {'Left', 'Right'};
    %stim = {'rs','rs','rs','rs','rs','Cont','Cont','Cont','Cont','Cont','Cont','Cont','Cont','Cont','Cont','Cont', '2Hz','2Hz','2Hz','2Hz','2Hz','2Hz','2Hz','2Hz','2Hz','2Hz'};
    %side = {'rs','rs','rs','rs','rs','Right','Right','Right','Right','Right','Right','Left','Left','Left','Left','Left', 'Right','Right','Right','Right','Right','Left','Left','Left','Left','Left'};

    for m = 1:length(mice)
        mouseID = mice{m};
        
%         allFiles = [];
%         for b = 1:length(BPOD_Files)
%             filePathBPOD = ['C:\Users\lmeyerb\OneDrive - Emory University\BPOD\',mouseID,'\Imaging_Whisker_Stim_Sweep_Frequency\Session Data\'];
%             stimSide = side{b};
%             % load in all the BPOD data for this session to save along with the
%            % imaging data
%            fileList = getBPODData(mouseID, filePathBPOD, BPOD_Files{b}, stimSide);
%            allFiles = [allFiles , fileList];
%         end

        for t = 1:length(trials)
            trial = trials(t);
            if trial <= 5
                tNum = ['00', num2str(trial)]; 
            else
                tNum = num2str(trial);
            end
            %seq_info.temp_reg = 0;
            out = preProcVSFP10L_Final(mouseID, '0918', tNum, mouseID); %,seq_info);
            %stimType  = stim{t};
            if t <= 50
                stimSide = 'Left';
            else
                stimSide = 'Right';
            end
            FolderName = ['X:\labs\keilholz-lab\Lisa\VSFP ButterFly\Data\VSFP_WhiskerStim\', mouseID, '_0918_', tNum,'.mat'];
            save(FolderName, 'out','stimSide')
        end

    end

end


function fileList = getBPODData(mouseID, filePathBPOD, BPOD_File, stimSide)

    load([filePathBPOD, BPOD_File])
    % load the Bpod file and store the trial time stamps
    RawEvents = SessionData.RawEvents.Trial;
    
    TrialStartTimestamp = SessionData.TrialStartTimestamp;
    TrialEndTimestamp = SessionData.TrialEndTimestamp;
    %change format of the timestamps to mm:ss
    TrialStartTimestamp = seconds(TrialStartTimestamp);
    TrialStartTimestamp.Format = 'mm:ss';
    TrialEndTimestamp = seconds(TrialEndTimestamp);
    TrialEndTimestamp.Format = 'mm:ss';

    fileName = erase(BPOD_File, '.mat');
    fileStartT = str2double(string(fileName((end-5):end)));
    fileStarT_str = string(fileName((end-5):end));
    
    fileDate = str2double(string(fileName((end-14):(end-7))));
    
    %convert the string to datetime format for easier arithmetic

    fileStartT_dt = datetime(fileStarT_str, 'InputFormat', 'HHmmss');
    fileStartT_dt = datetime(fileStartT_dt, 'Format', 'HH:mm:ss');
    
    TrialStartTimestamp = fileStartT_dt + TrialStartTimestamp;
    TrialEndTimestamp = fileStartT_dt + TrialEndTimestamp;
   
    % convert datetime to numeric value
    TrialStartTimestamp = str2double(string(datetime(TrialStartTimestamp, 'Format', 'HHmmss')));
    TrialEndTimestamp = str2double(string(datetime(TrialEndTimestamp, 'Format', 'HHmmss')));
    %go through each trial to find trigger time from MiCAM for each trial
    %for trial 
    totalTrialN = length(SessionData.RawEvents.Trial);
    ImgSyncTriggerTimes = zeros(1, totalTrialN);
    
    for trialN = 1:totalTrialN
        % If the MiCAM imaging happend, it should send a trigger to
        % Bpod via BNC1
        if isfield(RawEvents{1, trialN}.Events, 'BNC1High')
            imgSynTriggerT = RawEvents{1, trialN}.Events.BNC1High;
            ImgSyncTriggerTimes(1, trialN) = imgSynTriggerT(1);
            %BNC1 receives two triggers, use the first one as the
            %beginning of imaging
        else
            ImgSyncTriggerTimes(1, trialN) = nan;
        end
    
    
        % if all trials have received TTL triggers from Bpod 
        % every trial should be matched with imaging
        % if not, something is wrong
    
        numNonMatchedBpodTrials = sum(isnan(ImgSyncTriggerTimes));

        fileList(trialN).mouseID = mouseID;
        fileList(trialN).trialNum = trialN;
        fileList(trialN).StimSide = stimSide;
        fileList(trialN).fileDate = fileDate;
        fileList(trialN).TrialStartTimestamp = TrialStartTimestamp(trialN);
        fileList(trialN).TrialEndTimestamp = TrialEndTimestamp(trialN);
        fileList(trialN).ImgSyncTriggerTimes = ImgSyncTriggerTimes(trialN);
       
    end
end