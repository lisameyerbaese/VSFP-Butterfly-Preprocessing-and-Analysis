% Function takes pupil video file names stored in VSFP_50Hz_Proc.xlsx and
% calculates the pupil diameter using pupil_processing.m and then saves the
% generated .mat file in VSFP ButterFly/Data/Pupil Diameter
%
% Written by Lisa Meyer-Baese

function pupil_vidProc()

   %get root of file name based on computer being used 
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end

    addpath([startFile,'\VSFP ButterFly\Info'])
    
    T = readtable('VSFP_50Hz_Proc.xlsx');
    videos = T.video_proc_file;
    [videos,ind] = unique(videos);
    mice = T.Mouse;
    dates = T.Date;
    trials = T.Trials;
    
    %loop through all 50Hz videos, process them and save to \Data\Pupil
    %Diameter
    for i = 14:15%size(videos)
        k = ind(i);
        mouse = mice{k,1}; 
        date = dates(k,1); 
        trial = trials(k,1);
        vid_name = videos{i,1};
        
        outData = pupil_processing([],1,1);
        saved_data = strcat([startFile, '\VSFP ButterFly\Data\Pupil Diameter\',vid_name]);
        %save(saved_data,'outData')
    end
end