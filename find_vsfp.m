function [fpath, fstr, pathName] = find_vsfp(fileName, mouseID)

%%%%%%%%%%%%%%%%%%
%
%   function to find vsfp data that may be saved in various directories
%   Called within functions:
%       -readCMOS6.m
%
%   May need to update paths list to find VSFP files...
%       *also important to note sequence of paths for VC*
%
%%%%%%%%%%%%%%%%%%

%% Begin Search For VSFP Files 

% Get date from filename 
%fDate = [mouseID,'_',fileName(8:11),'_','2023'];  % for my collected whisker stim data
fDate = [mouseID,'_',fileName(9:12),'_','2018']; % for rhett's data
%fDate = [mouseID,'_',fileName(4:7),'_','2018']; % for rhett's data VSFP 30 +


%fName = [mouseID ,'_0', fileName(10:12),'_2018'];
%fName = [mouseID ,'_0', fileName(5:7),'_2018'];
[~, name] = system('hostname');
if contains(name,'jaeger')
    startFile = 'X:\labs\keilholz-lab\Lisa';
else
    startFile = 'X:\keilholz-lab\Lisa';
end

% Potential Directories that may hold file of interest
paths = {
%     ['/Volumes/My Book/JaegerLab/' mouseID '_' fDate '_2016/'],...
%     ['/Volumes/MAC_ONLY/' mouseID '/'],...
%     ['/Volumes/SD_0016/DJ_Lab/vsfp_imaging/' mouseID '_' fDate '_2016/'],...
%     ['/Volumes/SD_0016/DJ_Lab/vsfp_imaging/' mouseID '_' fDate '_2017/'],...
%     ['/Volumes/SD_0016/DJ_Lab/vsfp_imaging/' mouseID '_' fDate '_2018/'],...
%     ['/Volumes/SD_0016/DJ_Lab/vsfp_imaging/' mouseID '_' fDate '_2017/A_B/'],...
%     ['~/OneDrive - Emory University/' mouseID '/'],...
%     ['~/OneDrive - Emory University/' mouseID '_' fDate '_2016/'],...
%     ['~/OneDrive - Emory University/' mouseID '_' fDate '_2017/' mouseID '/'],...
%     ['/Volumes/My Book/New Folder/M2_1200Deg_CD1_Left/'],...
%     ['/Volumes/My Book/New Folder/M2_FreqSweep_CD1_Left/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/' mouseID '_' fDate '_2017/' mouseID '/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/VSFP17_0222_2018/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/VSFP17_0329_2018/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/VSFP12_1102_2016/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/' mouseID '_' fDate '_2018/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/' mouseID '_' fDate '_2017/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/' mouseID '_' fDate '_2016/'],...
%     ['/Volumes/djlab1_root/RawData/vsfp_raw/VSFP_Imaging/' mouseID '_' fDate '_2018/' mouseID '/'],...
%     ['/Volumes/WD_0018/Data/Personal/imaging/raw_data/' mouseID '_' fDate '_2018/'],...
%     ['/Volumes/WD_0018/Data/Personal/imaging/raw_data/' mouseID '_' fDate '_2017/']};
%    ['C:\Users\ywan848\Desktop\Data_examples_Jaeger_Lab_Miao\Cortical_imaging_Data\VSFP28_0801_2018\']};
%    ['Y:\keilholz-lab\Lisa\VSFP ButterFly\DIMG2\'], 
%['/Volumes/RhettData/Rhett_Voltage_Imaging-2020/vsfp_raw/VSFP_Imaging/VSFP24_0801_2018/'],
%['X:\keilholz-lab\Lisa\VSFP ButterFly\Ephys\VSFP26_0503_2018\'],
%['X:\keilholz-lab\Lisa\VSFP ButterFly\VSFP_LFP\Imaging\VSFP_Ephys3\VSFPE3_1216_2021\'],
%['D:\Rhett_Voltage_Imaging-2020\vsfp_raw\VSFP_Imaging\VSFP24_0801_2018\']};
%[startFile, '\VSFP ButterFly\VSFP_LFP\Imaging\VSFP_Ephys3\VSFPE3_1216_2021\']};
%[startFile, '\VSFP ButterFly\VSFP_LFP\Imaging\VSFP_Ephys4\VSFPE4_1216_2021\'];
['E:\Rhett_Voltage_Imaging-2020\vsfp_raw\VSFP_Imaging\', fDate, '\']};
%['D:\Raw Data\Voltage Imaging\', mouseID,'\', fDate, '\']};
%[startFile, '\VSFP ButterFly\RhettData\VSFP24_0801_2018\']};
%[startFile,'\VSFP ButterFly\VSFP_LFP\Imaging\VSFP_Ephys3\VSFPE3_1216_2021\']};


fids = zeros(1,length(paths)); 
x = 1;
done = 0;

while x < length(paths)+1
    if done == 0
    try
        fpath = [paths{x} fileName]
% Attempt to read the file
        fids(x)=fopen(fpath,'r','b');    
        fstr=fread(fids(x),'int8=>char')';
        if fids(x) > 0
            done = 1;
            disp(['Found file in ' paths{x} ' directory'])
        end
    catch
        %disp(['Unable to find file in ' paths{x} ' directory']) 
        x = x + 1;
    end
    y = x;
    else
        break
    end
end

pathName = paths{y};
% Attempy to read the file again
fclose(fids(y));

