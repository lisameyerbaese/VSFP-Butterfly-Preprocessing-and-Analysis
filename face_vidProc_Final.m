% Function takes face video file names stored in VSFP_50Hz_Proc.xlsx and
% calculates the amount of orofacial movements using pupil_processing.m 
% and then saves the generated .mat file in VSFP ButterFly/Data/Face
% States
%
% Written by Lisa Meyer-Baese

function face_vidProc_Final()

[~, name] = system('hostname');
if contains(name,'jaeger')
    startFile = 'X:\labs\keilholz-lab\Lisa';
else
    startFile = 'X:\keilholz-lab\Lisa';
end

addpath([startFile, '\VSFP ButterFly\Info'])


T = readtable('VSFP_50Hz_Proc.xlsx');
videos = T.video_proc_file;
[videos,ind] = unique(videos);
mice = T.Mouse;
dates = T.Date;
trials = T.Trials;

%loop through all 50Hz videos, process them and save to \Data\Pupil
%Diameter
for i = 14: size(videos)
    k = ind(i);
    mouse = mice{k,1}; 
    date = dates(k,1); 
    trial = trials(k,1);
    vid_name = videos{i,1};
    
    faceData = face_processing(0,0);
    saved_data = strcat([startFile, '\VSFP ButterFly\Data\Face States\',vid_name]);
    save(saved_data,'faceData')
%     saved_data = strcat('X:\keilholz-lab\Lisa\VSFP ButterFly\Data\Face States\','test');
%     save(saved_data, 'face_state')
end
end