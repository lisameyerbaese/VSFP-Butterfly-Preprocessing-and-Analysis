% Function takes xlsx sheet and generates figures for following analysis
%   1) takes all the traces, pupil, orofacial, hemo and volt and plots it
%       as a movie
%
% Inputs: 
%        optional mice number, or list of numbers i.e.
%        VSFP_MakeMove_OrofacialMvmt('VSFP24','VSFP25') or VSFP_MakeMove_OrofacialMvmt('VSFP24')
%        if no input is given, function runs through all animals
%
% Outputs: 
%        movie files for each of the trials 
%
% Written by Lisa Meyer-Baese

function VSFP_MakeMove_OrofacialMvmt(varargin)

    
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
    
    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_MakeMove_OrofacialMvmt\'];   % Your destination folder

    % select image used for mask
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    maskImg = imresize(maskImg, [100,100]); 
    maskImg = maskImg(:,:,1);
    maskImg = maskImg > 200;
   
    for m = 1: length(all_mice) 
        mouse = all_mice{m};
        ind = find(contains(T.Mouse,mouse));
        all_trials = T.Trials(ind);
        all_dates = T.Date(ind);
        videoFolder  = [startFile, '\VSFP ButterFly\Data\vsfp_50Hz_raw_pupil\'];

        % load in avg corr map to select ROIs
        disp('Checking to see if .mat Avg Corr Map already exists')
        fname = [mouse,'_','avgCorrMap.mat'];
        dataFolder  = [startFile, '\VSFP ButterFly\Data\VSFP_CorrMap_Pupil'];
        if isfile(fullfile(dataFolder, fname))
           disp('File Exists!')
           load(fullfile(dataFolder, fname), 'avgTrial_Volt');
        else 
           disp('File does not exists, running VSFP_CorrMap_Pupil.m to generate needed file')
           VSFP_CorrMap_Pupil(mouse)
           load(fullfile(dataFolder, fname), 'avgTrial_Volt');
        end
        
        f1 = figure(1);
        avgTrial_Volt = mean(avgTrial_Volt,3);
        imagesc(avgTrial_Volt), title('Select 4 ROI'), colorbar
        ROI1 = drawcircle('Color','k','FaceAlpha',0.4,'Label','1');
        ROI2 = drawcircle('Color','k','FaceAlpha',0.4, 'Label','2');
        ROI3 = drawcircle('Color','k','FaceAlpha',0.4, 'Label','3');
        ROI4 = drawcircle('Color','k','FaceAlpha',0.4, 'Label','4');
        
        % save all masks
        ROIs = cat(3, createMask(ROI1), createMask(ROI2), createMask(ROI3), createMask(ROI4));
        
        for k = 1:1:length(all_trials)
            %get corresponding pupil data
            FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
            date = FindTable.Date;
            trial = FindTable.Trials;
            videoFile = [FindTable.video_proc_file{1}, '.avi'];
            framesInd = FindTable.start_frame:1:FindTable.end_frame;
            header = [mouse ' ' num2str(date) ' ' num2str(trial) ' ' num2str(k)];
            
            % load in the corresponding face video file
            if isfile(fullfile(videoFolder, videoFile))
               disp('Video Exists!')
               movieObj = VideoReader([videoFolder, videoFile]);
            else 
               disp('Video file does not exists, skipping trial')
               break
            end

            % load the imaging data
            image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
            image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
            vsfp_data=load(image_data);

            % get pupil data
            pupil_data2x = vsfp_data.pupil_data2x; 
            pupil_data2x(1:2) = pupil_data2x(3); % first value has large jump

            % get face data and match the lengths
            face_data2x = reshape(vsfp_data.face_data2x, 1,[]); 
            lenPupil = min(length(face_data2x),length(pupil_data2x));
            face_data2x = face_data2x(1:lenPupil);
            pupil_data2x = pupil_data2x(1:lenPupil);
            
            %voltage signal 
            dataVolt = vsfp_data.out.imgDR3(:,:,1:lenPupil);
            % hemo signal 
            dataHemo = -1 * vsfp_data.out.hemoLP(:,:,1:lenPupil);
            
            % load dim of data 
            sXY = vsfp_data.out.sX * vsfp_data.out.sY;
            timeVec = vsfp_data.out.time(1:lenPupil);

            %load in data for each of the ROIs
            for r = 1:4
                % apply mask
                tempVolt = imrotate(dataVolt, 270) .* ROIs(:,:,r);
                tempHemo = imrotate(dataHemo, 270) .* ROIs(:,:,r);
                % Reshape data temp to avg traces for ROI
                tempVolt2D = mean(reshape(tempVolt,[sXY,lenPupil]),1);
                tempHemo2D = mean(reshape(tempHemo,[sXY,lenPupil]),1);
                % save data for plotting
                voltROI(r,:) = tempVolt2D;
                hemoROI(r,:) = tempHemo2D;
            end
            
            % make plot of data
            f2 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,2,1)
            imagesc(vsfp_data.out.baseD)
            title('Video')
            axis square, axis off

            subplot(2,4,5)
            imagesc(imrotate(vsfp_data.out.imgDS(:,:,200), 270) .* maskImg)
            %imagesc(imrotate(vsfp_data.out.imgD(:,:,20)/maxTempD, 270))
            hold on 
            viscircles(ROI1.Center, ROI1.Radius)
            viscircles(ROI2.Center, ROI2.Radius, 'color',"#EDB120")
            viscircles(ROI3.Center, ROI3.Radius, 'color',"#0072BD")
            viscircles(ROI4.Center, ROI4.Radius, 'Color',"#77AC30")
            title('Volt Signal')
            axis square, axis off
            caxis([0, 1])

            subplot(2,4,6)
            imagesc(imrotate(vsfp_data.out.hemoLP(:,:,200), 270) .* maskImg)
            %maxTempA = max(max(max(vsfp_data.out.imgA)));
            %imagesc(imrotate(vsfp_data.out.imgA(:,:,20)/maxTempA, 270))
            hold on 
            viscircles(ROI1.Center, ROI1.Radius)
            viscircles(ROI2.Center, ROI2.Radius, 'color',"#EDB120")
            viscircles(ROI3.Center, ROI3.Radius, 'color',"#0072BD")
            viscircles(ROI4.Center, ROI4.Radius, 'Color',"#77AC30")
            title('Hemo Signal')
            axis square, axis off
            caxis([0, 1])

            % plot traces as a function of time
            subplot(6,2,2)
            plot(timeVec, pupil_data2x,'color',"#D95319")
            h = gca;
            h.YAxis.Visible = 'off';
            title('Pupil Diameter')

            subplot(6,2,4)
            plot(timeVec, face_data2x, 'Color', "#7E2F8E")
            h = gca;
            h.YAxis.Visible = 'off';
            title('Orofacial Movement')

            subplot(6,2,6)
            plot(timeVec,voltROI(1,:), 'k')
            hold on 
            h = gca;
            h.YAxis.Visible = 'off';
            plot(timeVec,hemoROI(1,:), 'Color',"#0072BD")
            title('ROI 1','color','r')

            subplot(6,2,8)
            plot(timeVec,voltROI(2,:), 'k')
            hold on 
            h = gca;
            h.YAxis.Visible = 'off';
            plot(timeVec,hemoROI(2,:), 'Color',"#0072BD")
            title('ROI 2', 'color',"#EDB120")

            subplot(6,2,10)
            plot(timeVec,voltROI(3,:), 'k')
            hold on 
            h = gca;
            h.YAxis.Visible = 'off';
            plot(timeVec,hemoROI(3,:), 'Color',"#0072BD")
            title('ROI 3', 'color',"#0072BD")

            subplot(6,2,12)
            plot(timeVec,voltROI(4,:), 'k')
            hold on 
            h = gca;
            h.YAxis.Visible = 'off';
            xlabel('Time (s)')
            plot(timeVec,hemoROI(4,:), 'Color',"#0072BD")
            title('ROI 4', 'Color',"#77AC30")
          
            fname = [mouse,'_0', num2str(date) ,'_',num2str(trial),'_TestFig.fig'];
            saveFig(f2, fname, FolderName);
            close (f2);

            % make a movie of the data
            displayMotion(movieObj,ROI1, ROI2,ROI3,ROI4,vsfp_data, pupil_data2x,face_data2x, timeVec,voltROI, hemoROI,framesInd, header, maskImg)
        end
    end   
end

function displayMotion(movieObj,ROI1, ROI2,ROI3,ROI4,vsfp_data, pupil_data2x,face_data2x, timeVec,voltROI, hemoROI,framesInd, header, maskImg)

    file = strcat(header, '.mp4');
    vidfile = VideoWriter(file,'MPEG-4');
    open(vidfile);
    
    ind = 1;
    % update the video file
    for frame = 1:2:length(pupil_data2x)
        f3 = figure('units','normalized','outerposition',[0 0 1 1]);
        %set(f3, 'Visible', 'off');
        % make plot of data, ind is used to account for the fact that video
        % data is sampled at a slower rate
        tempFrame = framesInd(ind);
        ind = ind + 1;
        thisFrame = read(movieObj, tempFrame);
        %subplot(2,3,1)
        subplot(2,4,1)
        image(thisFrame)
        set(gca,'ytick',[], 'xtick',[])
        axis square
        title('Video')
        axis square, axis off

       subplot(2,4,2)
       rawData = imrotate(vsfp_data.out.imgD(:,:,frame) .* 10, 270) .* maskImg;
       imagesc(rawData, [0, 100000])
       title('Donor Channel Raw Signal')
       axis square, axis off

       subplot(2,4,5)
       voltData = imrotate(vsfp_data.out.imgDR3(:,:,frame) .* 10, 270) .* maskImg;
       %dData = imrotate(vsfp_data.out.imgD(:,:,frame), 270) .* maskImg;
       imagesc(imgaussfilt(voltData,0.5))
       hold on 
       viscircles(ROI1.Center, ROI1.Radius)
       viscircles(ROI2.Center, ROI2.Radius, 'color',"#EDB120")
       viscircles(ROI3.Center, ROI3.Radius, 'color',"#0072BD")
       viscircles(ROI4.Center, ROI4.Radius, 'Color',"#77AC30")
       title('Voltage Signal')
       axis square, axis off
       clim([-0.2, 0.2])
       %clim([0, 5000])
      
       subplot(2,4,6)
       hemoData = imrotate(-1* vsfp_data.out.hemoLP(:,:,frame) .* 10, 270) .* maskImg;
       %aData = imrotate(vsfp_data.out.imgA(:,:,frame), 270) .* maskImg;
       imagesc(imgaussfilt(hemoData,0.5))
       hold on 
       viscircles(ROI1.Center, ROI1.Radius)
       viscircles(ROI2.Center, ROI2.Radius, 'color',"#EDB120")
       viscircles(ROI3.Center, ROI3.Radius, 'color',"#0072BD")
       viscircles(ROI4.Center, ROI4.Radius, 'Color',"#77AC30")
       title('Hemo Signal')
       axis square, axis off
       clim([-0.2, 0.2])
       %clim([0, 5000])

       % plot traces as a function of time
       subplot(6,2,2)
       plot(timeVec, pupil_data2x,'color',"#D95319");
       hold on
       h = gca;
       h.YAxis.Visible = 'off';
       xline(timeVec(frame),'-.')
       title('Pupil Diameter')

       subplot(6,2,4)
       plot(timeVec, face_data2x, 'Color', "#7E2F8E")
       hold on 
       xline(timeVec(frame),'-.')
       h = gca;
       h.YAxis.Visible = 'off';
       title('Orofacial Movement')

       subplot(6,2,6)
       plot(timeVec,voltROI(1,:), 'k')
       hold on 
       h = gca;
       h.YAxis.Visible = 'off';
       xline(timeVec(frame),'-.')
       plot(timeVec,hemoROI(1,:), 'Color',"#0072BD")
       title('ROI 1','color','r')

       subplot(6,2,8)
       plot(timeVec,voltROI(2,:), 'k')
       hold on 
       h = gca;
       h.YAxis.Visible = 'off';
       xline(timeVec(frame),'-.')
       plot(timeVec,hemoROI(2,:), 'Color',"#0072BD")
       title('ROI 2', 'color',"#EDB120")

       subplot(6,2,10)
       plot(timeVec,voltROI(3,:), 'k')
       hold on 
       h = gca;
       h.YAxis.Visible = 'off';
       xline(timeVec(frame),'-.')
       plot(timeVec,hemoROI(3,:), 'Color',"#0072BD")
       title('ROI 3', 'color',"#0072BD")

       subplot(6,2,12)
       plot(timeVec,voltROI(4,:), 'k')
       hold on 
       h = gca;
       h.YAxis.Visible = 'off';
       xlabel('Time (s)')
       xline(timeVec(frame),'-.')
       plot(timeVec,hemoROI(4,:), 'Color',"#0072BD")
       title('ROI 4', 'Color',"#77AC30")

       sgtitle(header)
        
       %add to video file
       drawnow; 
       F(frame) = getframe(gcf); 
       writeVideo(vidfile,F(frame));
       pause(0.1)
       close(f3)
    end
    close(vidfile)
end