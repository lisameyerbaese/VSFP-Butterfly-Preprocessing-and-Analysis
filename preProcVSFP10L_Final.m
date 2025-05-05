function out = preProcVSFP10L_Final(fpre, fDate, fNum, mouseID, seq_info) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SWITCHED CAMERA A AND CAMERA B TO BE FOR ACCEPTOR AND
% DONAR CHANNELS RESPECTIVLEY
%   
%
% Preprocessing for VSFP imaging data (Dual Camera Acquisition) 
% 
% Usage:
%   out = preProcVSFP(fDate, fNum, avgGain)
% where: avgGain = 0 or 1 (calculate gain factors from trial avg)
%
% 7/16/2015 - Much faster Version 2.0
% 9/02/2015 - Included "method" for specifiying method of removing HR
% artifact & "stimTrial" specification for removing stim artifact when
% performing PCA
% 10/26/2015 - Script overhaul and simplification, includes option to 
% generate equalization gain factor from trial averages (per conv with T. Knopfel)
% 12/09/2015 - Added compatibility for DIF acquired images through ultima
% dual camera system (critical for ratiometric equalization)
% 05/13/2016 - When SAVE is turned on (by setting save equal to 1 in bottom
% of the script), preprocessed output will be saved as 'VSFP_Output_
% 07/18/2016 - Modified stim removal function to rescale stimulus shifted signal
% by comparing STD of Hemodynamic signals pre and during stimulation
% 08/05/2018 - Updated spatialAvg function (spatialAvg2) that uses conv2 instead of
% slower matrix averaging (updated from version 7)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test program with test data set *** From previous version ***
% if strcmp(fDate,'test')
%     fDate = 'test';
%     fNum = '000';
%     mouseID = 'test';
%     [imgA, imgD] = getTestData;
%     imgA1 = imgA;
%     imgD1 = imgD;
% else

%% Save File ON - Not currently in use
% saveVar = 1;

%% Load sequence info
if exist('seq_info','var')
    disp('Loaded info from previous sequence...')
    % Check for filename info
    if isfield(seq_info,'fname')
        disp(['Info from file ' seq_info.fname ' ...'])
    end
    % Check for ref_imgs 
    if isfield(seq_info,'ref_imgs')
        ref_imgs = seq_info.ref_imgs;
    end
    % Check for avgGain 
    if isfield(seq_info,'avgGain')
        avgGain = seq_info.avgGain;
    end
    % Check for 
    if isfield(seq_info,'ref_masks')
        ref_masks = seq_info.ref_masks;
    end
    if isfield(seq_info,'temp_reg')
        temp_reg = 1;
    else
        temp_reg = 0;
    end
    if isfield(seq_info,'norm_method')
        norm_method = seq_info.norm_method;
    else
        norm_method = 'trial';
    end
end

%% Load imaging files 
% file prefix
% Load the prefix from the a header file in the directory

%fpre = 'VSFP_01A0'; % for Rhetts original data
%fpre = 'VSFPE4_';
%fpre = 'VSFP_01A0';
%fpre = mouseID;  %for my whisker stim data

% Convert inputs if needed
fStr = num2str(fNum);

if length(fStr) == 2
    fStr = ['0' fStr];
elseif length(fStr) < 2
    fStr = ['00' fStr];
end 
fNum = str2double(fStr);
disp(['Processing file ' (fStr) ' ...'])

% Load image filesExcessData, pathName, StmTrg, fpath
%[imgA, anaA1, anaA2, Astim, pathName, StmTrg, fpath,~,Fs] = readCMOS6([fpre,'_', num2str(fDate) '-' fStr '_A.rsh'], mouseID); % for whisker stim and ephys data
[imgA, anaA1, anaA2, Astim, pathName, StmTrg, fpath,~,Fs] = readCMOS6([fpre, num2str(fDate) '-' fStr '_A.rsh'], mouseID); % for rhett's data
%[imgA, anaA1, anaA2, Astim, pathName, StmTrg, fpath,~,Fs] = readCMOS6(['1000', num2str(fDate) '-' fStr '_A.rsh'], mouseID); % for rhett's data VSFP30 +
%[imgD] = readCMOS6([fpre,'_', num2str(fDate) '-' fStr '_B.rsh'],mouseID); % for whisker stim and Ephys data
[imgD] = readCMOS6([fpre, num2str(fDate) '-' fStr '_B.rsh'],mouseID); % for rhett's data
%[imgD] = readCMOS6(['1000', num2str(fDate) '-' fStr '_B.rsh'],mouseID); % for rhett's data VSFP30 + 

% Values are currently assigned by stimCheck function
% imgA1 = imgA;
% imgD1 = imgD;

%% Check image acquisition method (CDF vs DIF)

if imgA(50,50,1) < 1000
    DIF = 1;
    disp('DIF used to acquire images! Recalculating with base fluoresence...')

% for image A   
    %pathA = [pathName '1000' num2str(fDate) '-' fStr 'A.rsm'];
    %pathA = [pathName fpre '_' num2str(fDate) '-' fStr 'A.rsm']; % whisker stim data
    pathA = [pathName fpre ,num2str(fDate) '-' fStr 'A.rsm']; % Rhetts data
    %pathA = [pathName '1000', num2str(fDate) '-' fStr 'A.rsm']; % Rhetts data VSFP 30+
    fidA = fopen(pathA,'r','n');
    fdataA = fread(fidA,'int16');
    fdataA = reshape(fdataA,128,100);
    baseA = fdataA(21:120,:)';
    imgA = bsxfun(@plus,imgA,baseA);
    baseAre = reshape(baseA,[100,100]);

% for image D    
    %pathD = [pathName '1000' num2str(fDate) '-' fStr 'B.rsm'];
    %pathD = [pathName fpre '_' num2str(fDate) '-' fStr 'B.rsm']; % whisker stim data
    pathD = [pathName fpre ,num2str(fDate) '-' fStr 'B.rsm']; % Rhett's data
    %pathD = [pathName '1000', num2str(fDate) '-' fStr 'B.rsm']; % Rhetts data VSFP 30+
    fidD = fopen(pathD,'r','n');
    fdataD = fread(fidD,'int16');
    fdataD = reshape(fdataD,128,100);
    baseD = fdataD(21:120,:)';
    imgD = bsxfun(@plus,imgD,baseD);
    baseDre = reshape(baseD,[100,100]);
    disp('done!')
else
    DIF = 0;
end

%% Register images (Now 2 parts of registration)
% 1. template image for each mouse (rigid registration rotational + translational) 
% 2. inter-trial registration (subpixel translational)
if exist('seq_info','var')
    if temp_reg == 1 % load template image for corresponding mouse
        disp('Template Registration...')
        addpath( 'X:\labs\keilholz-lab\Lisa\VSFP ButterFly\Data\VSFP_WhiskerStim_ImagingData\')
        temp_load = load([mouseID '_reg.mat']);
        temp_reg = temp_load.reg;   
        temp_base = temp_reg.baseImg;
    
        [optimizer, metric] = imregconfig('monomodal'); % monomodal when all cameras are donor channel
%     baseD_tempReg = imregister(baseD, temp_base, 'rigid', optimizer, metric);
        tform_tempReg = imregtform(baseD,temp_base,'rigid',optimizer,metric);
    
        [imgA,imgD] = tempreg_vsfp(imgA,imgD,tform_tempReg);
        baseD_tempReg = imgD(:,:,1);
        baseA_tempReg = imgA(:,:,1);
        ref_imgs(:,:,1) = baseA_tempReg;
        ref_imgs(:,:,2) = baseD_tempReg;
        disp('done!')
    end
end

%% inter-trial registration

if 1
    disp('Registering images...')
    if exist('ref_imgs','var') == 1 
        if ~isempty(ref_imgs)
            [imgA,imgD,reg_info] = dftreg_vsfp(imgA,imgD,100,ref_imgs(:,:,1),ref_imgs(:,:,2));
            imgReg = 1;
        end
        imgReg = 0;
    else
        [imgA,imgD,reg_info] = dftreg_vsfp(imgA,imgD,100,baseA,baseD);
        imgReg = 1;
    end
end

%% Create Mask of Pixels (based on fluorescent intensity)

% mask = zeros(10000,1);
% aRe = reshape(imgA(:,:,1),[10000,1]);
% mask((aRe > 1500)) = 1;
% mask = reshape(mask,[100,100]);

% Dimensions of the loaded files?
if isequal(size(imgD),size(imgA)) == 1
    [sX, sY, sZ] = size(imgD);
    sXY = sX*sY;
else
    error('Donor and Acceptor Files must be the same size')
end

% How many frames are in each sequence (for me usually 2048)
numFrames = sZ;

% Time array (sampling at 5 msec = 0.005 sec)
FsImg = Fs; % Sampling Rate (Hz)
disp(['Sampling freq: ' num2str(FsImg) ' ...'])
time = (1:numFrames).*(1/FsImg); % array in seconds

% Select ROI encompassing entire brain FOV and create masks
% !This is old from version before but may get added back in future
% version...
% if exist('region','var') ~= 1 && exist('PCA','var') == 1
%     region = roiSelect(fDate, fStr);
% else
%     region = [1 100 1 100];
% end

Dphys = imgD;
Aphys = imgA;

% Number of frames to use for baseline average -- typically first 200 frames (1sec) of data
if exist('norm_method','var')
    if strcmp(norm_method,'trial')
        avgingFr = numFrames;
    elseif strcmp(norm_method,'start') % normalize to first 500ms of activity
        avgingFr = floor(FsImg./5); % 100 frames at 200Hz, 25 frames at 50Hz
    else
        disp('did not recognize norm_method... using first 4 frames')
        avgingFr = 4; % normalize to first 4 frames - not recommended...
    end
else
    disp('norm_method not specified using seq_info... using first 200ms avg')
    avgingFr = floor(FsImg./5);
end

% Create stucture for output and version control
out = struct('fileNum',[num2str(fDate) '-' num2str(fNum)],'version','v3.0','date',date);

%% Equalized Ratio to correct for wavelength-dependent loght absorption by hemoglobin
% Calculate standard deviations of Acceptor and Donor channels 

% Can modify filter settings for optimal HR subtraction
if exist('avgGain','var')
    HRwin = [0.2*FsImg,0.25*FsImg];
else
    %HRwin = HRfilter2(imgD,fDate,fStr,FsImg,mouseID);
    HRwin = HRfilter3(imgD,fDate,fStr,FsImg,mouseID);
end
disp(['Hemodynamic range at ' num2str(HRwin(1)) ' and ' num2str(HRwin(2)) ' ...'])

% if exist('avgGain','var') == 1    
%     HRwin = [0.5 1.5];
% else
%     HRwin = [12 15];
% end

HRpeak = HRwin(2)-((HRwin(2)-HRwin(1))/2);
disp(['Hemodynamic peak at ' num2str(HRpeak) ' ...'])
disp(['Filtering between ' num2str(HRwin(1)) ' and ' num2str(HRwin(2)) ' ...'])
hemfilter = make_ChebII_filter(1, FsImg, HRwin, [HRwin(1)*0.9 HRwin(2)*1.1], 20);

% Preallocate couple variables for speed :)
Ahem2d = ones(sXY,numFrames);
Dhem2d = ones(sXY,numFrames);
blur3 = [];

%% First perform operations outside of cycle:
        
% Determine if presence of optogenetic stimulation
% [imgA1,imgD1,I] = stimCheck(imgD,Aphys,Dphys,sX,sY);
imgA1 = Aphys;
imgD1 = Dphys;
I = 0;
Dphys2d = reshape(imgD1, [sXY, sZ]);
Aphys2d = reshape(imgA1, [sXY, sZ]);
stdA2 = std(Aphys2d(:,:),0,2);
stdD2 = std(Dphys2d(:,:),0,2);  

if exist('ref_masks','var')
    maskA = squeeze(ref_masks(:,:,1));
    maskD = squeeze(ref_masks(:,:,2));
else
    maskA = imgA(:,:,1)>quantile(reshape(imgA(:,:,1),[1,10000]),0.50);
    maskD = imgD(:,:,1)>quantile(reshape(imgD(:,:,1),[1,10000]),0.50);
%     [maskA,maskD] = slider_mask(imgA(:,:,1),imgD(:,:,1));
end

% Averages of baseline VSFP recordings
Aavg = mean(Aphys2d(1:sXY,1:avgingFr),2);
Davg = mean(Dphys2d(1:sXY,1:avgingFr),2);

%% Cycle through each pixel (100 x 100 window) *This is slow! - remove in future versions*
for xy = 1:sXY
% filter for hemodynamic signal (12-14 hz) only for the length of the
% averaging frames ***
    Ahem2d(xy,:) = filtfilt(hemfilter.numer, hemfilter.denom,Aphys2d(xy,1:numFrames));
    Dhem2d(xy,:) = filtfilt(hemfilter.numer, hemfilter.denom,Dphys2d(xy,1:numFrames));
end 


% Standard deviation of filtered signal   
stdA = std(Ahem2d(:,:),0,2);
stdD = std(Dhem2d(:,:),0,2);  

% avgGain = avgGain

% gamma and delta values for each pixel    
if exist('avgGain','var') == 1
    disp('Using avg gain to calculate gamma/delta...')
%     gamma1 = avgGain(:,:,1).*physMask(:,:,1);
%     gamma = reshape(gamma1, [sXY, 1]);
      gamma = avgGain(:,1);
%     delta1 = avgGain(:,:,2).*physMask(:,:,1);
%     delta = reshape(delta1, [sXY, 1]);
      delta = avgGain(:,2);

else
    gamma = 0.5*(1+((Aavg.*stdD)./(Davg.*stdA)));
    delta = 0.5*(1+((Davg.*stdA)./(Aavg.*stdD)));
end

% Calculate equalized A and D data
imgDiffA = bsxfun(@minus,Aphys2d,Aavg);
imgDiffD = bsxfun(@minus,Dphys2d,Davg);

Aeql = bsxfun(@plus,bsxfun(@times,gamma,imgDiffA),Aavg);
Deql = bsxfun(@plus,bsxfun(@times,delta,imgDiffD),Davg);

% Calculate new baseline values for equalized frames
avgAe = mean(Aeql(:,1:avgingFr),2);
avgDe = mean(Deql(:,1:avgingFr),2);        
        
% and gain-corrected rationmetric signal
imgDiv = Aeql ./ Deql;
avgDiv = avgDe ./ avgAe;          
imgDR = bsxfun(@times,imgDiv,avgDiv)-1;
imgSum = Aeql .* avgDe + Deql .* avgAe;
imgDS = (imgSum-mean(imgSum(:,1:avgingFr),2))./mean(imgSum(:,1:avgingFr),2);

%% Regress out remaining hemodynamic activity from voltage signal
% mask hemodynamic signals
disp('regress out low-freq hemodynamic activity ...') 
mask2d = reshape(maskD,[1,sX*sY]);
imgDS_masked = imgDS(mask2d,:);

% zscore signals so that large artifacts do not dominate the average signal
imgDS_zscored = zscore(imgDS_masked')';
mean_imgDS = mean(imgDS_zscored);

% filter mean-masked hemodynamic activity between 1-5Hz
hem_lpfilt = make_ChebII_filter(2, FsImg, 5, 5*1.1, 20);
mean_imgDS_filt = filtfilt(hem_lpfilt.numer, hem_lpfilt.denom, mean_imgDS);
   
% Run gsr (regress global hemodynamic signal from each pixel between 1-5Hz)
imgDR3 = vsfp3_gsr(imgDR,mean_imgDS_filt);

% Mean of voltage signal
imgDR3_masked = imgDR3(:,mask2d);
imgDR3_zscored = zscore(imgDR3_masked);
mean_imgDR3 = mean(imgDR3_zscored');

% and transpose imgDR3
imgDR3 = imgDR3';

%
%% Spatial Blur VSFP data - Currently removed in Version 5!

% Blur now turned off for data - to be processed later
% blur = 0;
% if exist('blur','var') == 1
%     disp('Spatial Blurring of VSFP data... ')
%     blur3 = spatialAvg2(reshape(imgDR3,[sX,sY,sZ]),3);
%     disp('Spatial Blurring of Hem data... ')
%     hem3 = spatialAvg2(reshape(imgDS,[sX,sY,sZ]),3);
% end

%% Output Variables and Save
out.imgA = imgA;
out.imgD = imgD;

out.mean_imgDS = mean_imgDS;
out.mean_imgDS_filt = mean_imgDS_filt;
out.mean_imgDR3 = mean_imgDR3;
% out.imgD = imgD;
% out.imgA = imgA;
% out.imgA1 = imgA1;
% out.imgD1 = imgD1;
out.fNum = fStr;
out.fDate = fDate;
% out.imgDR2 = reshape(imgDR,[sX,sY,sZ]);
out.imgDR3 = reshape(imgDR3,[sX,sY,sZ]);
out.imgDS = reshape(imgDS,[sX,sY,sZ]);
out.time = time;
out.gamma = gamma;
out.delta = delta;
% out.imgVoltFilt = reshape(imgVoltFilt,[sX,sY,sZ]);
% out.imgDRFilt = reshape(imgDRFilt,[sX,sY,sZ]);
out.anaA1 = anaA1;
out.anaA2 = anaA2;
% out.stim1 = Astim;
% out.DIF = DIF;
out.baseA = baseAre;
out.baseD = baseDre;
% out.blur3 = blur3;
% out.hem3 = hem3;
% out.stdA = reshape(stdA,[sX,sY]);
% out.stdD = reshape(stdD,[sX,sY]);
% out.stdA2 = reshape(stdA2,[sX,sY]);
% out.stdD2 = reshape(stdD2,[sX,sY]);
% out.Aeql = reshape(Aeql,[sX,sY,sZ]);
% out.Deql = reshape(Deql,[sX,sY,sZ]);
% out.Aavg = reshape(Aavg,[sX,sY]);
% out.Davg = reshape(Davg,[sX,sY]);
% out.Ahem = reshape(Ahem2d,[sX,sY,sZ]);
% out.Dhem = reshape(Dhem2d,[sX,sY,sZ]);
% out.avgAe = reshape(avgAe,[sX,sY]);
% out.avgDe = reshape(avgDe,[sX,sY]);
% out.imgA1 = reshape(imgA1,[sX,sY,sZ]);
% out.imgD1 = reshape(imgD1,[sX,sY,sZ]);
% out.imgGD = reshape(imgGD,[sX,sY,sZ]);
% out.imgGA = reshape(imgGA,[sX,sY,sZ]);
% out.imgGAFilt = reshape(imgGDFilt,[sX,sY,sZ]);
% out.imgGDFilt = reshape(imgGAFilt,[sX,sY,sZ]);
out.version = 'PreProcVSFP10L_win';
out.StmTrg = StmTrg;
out.mask = maskD;
out.maskA = maskA;
out.maskD = maskD;
out.sX = sX;
out.sY = sY;
out.sZ = sZ;
% out.imgDiffA = reshape(imgDiffA,[sX,sY,sZ]);
% out.imgDiffD = reshape(imgDiffD,[sX,sY,sZ]);
out.mouseID = mouseID;
out.Fs = Fs;
out.pathName = pathName;
out.imgReg = imgReg;
if imgReg == 1
    out.regInfo = reg_info;
else
    out.regInfo = [];
end
out.HRwin = HRwin;

% if saveVar == 1
%     save(['VSFP_Out_' fStr '_' num2str(fDate) '.mat'],'out');
% end

end

