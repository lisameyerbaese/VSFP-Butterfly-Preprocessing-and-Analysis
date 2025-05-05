% Function takes xlsx sheet and generates figures for following analysis
%   1) loads in all the imaging data, and looks at the number of trials and
%       length of trials for each mouse 
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

function VSFP_DataDistribution(varargin)
    
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

    m1 = [];
    m2 = [];
    m3 = [];
    m4 = [];
    m5 = [];

    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        
        trialLen = [];
        for k = 1:length(all_trials)

            %get corresponding pupil data
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;

            % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            gsVolt = vsfp_data.out.gsVolt;
            Fs = vsfp_data.out.Fs;

            trialLen = [trialLen ,length(gsVolt) / Fs];
            
        end

        switch m 
            case 1
                 m1 = trialLen;
            case 2
                 m2 = trialLen;
            case 3
                 m3 = trialLen;
            case 4
                 m4 = trialLen;
            case 5
                 m5 = trialLen;
        end
    end

    display('Trial duration for each of the 5 mice is saved in M1 - M5 ')

end