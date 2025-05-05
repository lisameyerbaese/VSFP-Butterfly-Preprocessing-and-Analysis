function face_state = face_processing(plot_fig, plot_vid)
%%%
%
%   Processing Pipeline for Face Movement Detection: The onset and offset of movement periods were 
%   detected based on average changes in pixel intensity between frames. Mean orofacial movement 
%   signals were differentiated and squared, then thresholded. 
%   
%   Inputs: plot_fig = 1 will output figures ( = 0 will not generate
%                           plots)
%           plot_vid = 1  will output final video showing change in motion
%                          over time ( = 0 will not create video file)
%           Example face_state = face_processing(1, 1)
%   Outputs: face_states a structure will all the calculated data
%   Copyright MIT (c) 2021 Lisa Meyer-Baese
%%%

%% Loading the camera files
    close all;
    
    persistent lastPath pathname
    % If this is the first time running the function this session,
    % Initialize lastPath to 0
    if isempty(lastPath) 
        lastPath = 0;
    end

    if lastPath == 0 % First time calling 'uigetfile', use the pwd
        [filename, pathname] = uigetfile('*.avi', 'Select Video File');
    else % All subsequent calls, use the path to the last selected file
        [filename, pathname] = uigetfile([lastPath '/*.avi'], 'Select Video File');
    end
    % Use the path to the last selected file
    % If 'uigetfile' is called, but no item is selected, 'lastPath' is not overwritten with 0
    if pathname ~= 0
        lastPath = pathname;
    end

    if filename == 0
        error('No file selected')
    end
    disp(['-     loading video file: ' filename])
    if exist([pathname filename(1:end-4) '.mat'], 'file')
        user_in = input('File Exists - hit enter to continue: ');
        if user_in ~= 1
            error('Processing ended... enter "1" to re-analyse')
        end
    end
    movieObj = VideoReader([pathname, filename]);

    % Determine how many frames there are.
	numberOfFrames = movieObj.NumberOfFrames;
	vidHeight = movieObj.Height;
	vidWidth = movieObj.Width;
    fs = movieObj.FrameRate;
	
    faceData = zeros(vidHeight,vidWidth,numberOfFrames);
	numberOfFramesWritten = 0;

    %% manually crop video file to contain only nose, snout, jaw and whiskers 
    figure(1)
    image(read(movieObj, 1))  % read in the first frame to crop out all unwanted areas
    roi = drawrectangle('Label','ROI','Color',[1 0 0]);
    flag = 0;
    while flag == 0
        user_in = input('ROI selected: hit 1 to continue, hit 2 to try again ');
        if user_in == 1
            flag = 1; 
        elseif user_in ~= 1
           figure(1)
           image(read(movieObj, 1))  % read in the first frame to crop out all unwanted areas
           roi = drawrectangle('Label','ROI','Color',[1 0 0]);
        end
    end
    mask = createMask(roi);   % create a mask for region of interest, to apply and crop all other frames
    
    %%
	% Loop through the movie, writing all frames out.
	% Each frame will be in a separate file with unique name.
	meanGrayLevels = zeros(numberOfFrames, 1);

    for frame = 1 : numberOfFrames
		% Extract the frame from the movie structure.
		thisFrame = read(movieObj, frame);

		% Calculate the mean gray level.
		grayImage = double(rgb2gray(thisFrame));
        thisFrame = grayImage .* mask;  %apply mask to frame      
        faceData(:,:,frame) = thisFrame; %grayImage;
        
        progressIndication = sprintf('Processed frame %4d of %d.', frame, numberOfFrames);
        disp(progressIndication)   
    end

    og_faceData = faceData;
    
    %% face_data caclulated 
    faceData = squeeze(mean(mean(og_faceData)));
    times = 1/fs:1/fs:length(faceData).*1/fs;
    mvmt = diff(faceData).^2;   % calculates differences between adjacent elements of X along the first array dimension whose size does not equal 1
    mvmt_filt = medfilt1(mvmt,floor(fs),'truncate');
    mvmt_filt_z = zscore(mvmt_filt);
    mvmt_thr = 0; %1.5.*std(mvmt_filt_z);
    mvmt_inds = find(mvmt_filt_z>mvmt_thr);
    rest_inds = find(mvmt_filt_z<mvmt_thr);

    %% Find continuous periods of rest and movement
    mvmt_inds01 = zeros(1,length(mvmt_inds));mvmt_inds01 = zeros(1,length(mvmt_inds));
    mvmt_inds01(mvmt_inds) = 1;
    mvmt_inds01(1) = 0; mvmt_inds01(end) = 0;
    mvmt_pulses = findPulses(mvmt_inds01);

    rest_inds01 = zeros(1,length(rest_inds));
    rest_inds01(rest_inds) = 1;
    rest_inds01(1) = 0; rest_inds01(end) = 0;
    rest_pulses = findPulses(rest_inds01); 

    mvmt_periods = mvmt_pulses.ends-mvmt_pulses.starts;
    rest_periods = rest_pulses.ends-rest_pulses.starts;

    %% Histogram of rest vs movement period lengths
    edges = 1:5:500;
    if plot_fig == 1
        figure() 
        histogram(mvmt_periods,edges), hold on
        histogram(rest_periods,edges)
        title(' Histogram of Rest vs Movement Period Lengths')
        legend ('Movement', 'Rest')
    end
    %% Find periods greater than some time threshold (using 2s)
    min_sec = 1;
    max_sec = 4;
    min_length = min_sec.*fs;
    max_length = max_sec.*fs;
    mvmt_good = find(mvmt_periods>min_length);
    rest_good = find(rest_periods>min_length);
    %% Sample plot
    if plot_fig == 1
        figure() 
        plot(times(1:end-1), mvmt_filt_z)
        hold on
        line([times(1),times(end-1)],[mvmt_thr,mvmt_thr])
        scatter(times(mvmt_inds),ones(length(mvmt_inds),1).*max(mvmt_filt))
        title('Time Course with thresholded movements')
    end

    %% 
    mvmt_good_starts = mvmt_pulses.starts(mvmt_good);
    rest_good_starts = rest_pulses.starts(rest_good);
    
    % Check to make sure none of the ranges fall outside the length of data
    rest_good_starts((rest_good_starts-max_length/2)<1) = [];
    rest_good_starts((rest_good_starts+max_length/2)>length(mvmt)) = [];
    mvmt_good_starts((mvmt_good_starts-max_length/2)<1) = [];
    mvmt_good_starts((mvmt_good_starts+max_length/2)>length(mvmt)) = [];

    mvmt_onset_cnst_range = [mvmt_good_starts-max_length./2;mvmt_good_starts+max_length./2];
    rest_onset_cnst_range = [rest_good_starts-max_length./2;rest_good_starts+max_length./2];

    %% structure with all the good stuff
    face_state = struct();
    face_state.mvmt_cnst_range = mvmt_onset_cnst_range;
    face_state.mvmt_pulses = mvmt_pulses;
    face_state.rest_cnst_range = rest_onset_cnst_range;
    face_state.rest_pulses = rest_pulses;
    face_state.min_sec = min_sec;
    face_state.max_sec = max_sec;
    face_state.fs = fs;
    face_state.mvmt_filt_z = mvmt_filt_z;
    face_state.mvmt_z = zscore(mvmt); 
    
    
    %% Generate final video of results
     if plot_vid == 1
         header = strcat([filename,'face_processing']);  % name of video file
         time = (1:length(mvmt_filt_z))/fs;
         displayMotion(movieObj, mvmt_filt_z, header, time, fs)
     end
    
end


function displayMotion(movieObj, faceData, header, time, Fs)
    ds = 10;
    file = strcat(header, '.mp4');
    vidfile = VideoWriter(file,'MPEG-4');
    open(vidfile);
    
    %firstFrame = frames(1);
    
    figure()
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %plot data for the first frame
    subplot(2,2,[1 2])
    image(read(movieObj, 1))
    title(header)
    
    subplot(2,2,3)
    plot(time, faceData)
    hold on
    h1  = plot(1, faceData(1),'mo','MarkerFaceColor','m','YDataSource','Y','XDataSource','X');
    
    subplot(2,2,4)
    range = [1:1/Fs:10];
    plot(time,faceData)
    hold on 
    h2 = plot(1, faceData(1),'mo','MarkerFaceColor','m','YDataSource','Y','XDataSource','X');
    
    for frame = 1:ds : length(faceData)
        thisFrame = read(movieObj, frame);
        subplot(2,2,[1 2])
        image(thisFrame)
        set(gca,'ytick',[], 'xtick',[])
        axis square
        
        subplot(2,2,3)
        X = time(frame);
        Y = faceData(frame);
        refreshdata(h1,'caller');
        set(gca,'ytick',[])
        axis( [0 max(time) -2 (max(faceData) + 3) ])
        title('Change in Orofacial Movements Over Time')
        xlabel('Time (s)')
        
        subplot(2,2,4)
        X = time(frame);
        Y = faceData(frame);
        range = range + 1/Fs;
        axis([range(1), range(end) , faceData(frame) - 5 ,faceData(frame) + 5])
        %refreshdata(h2,'caller');
        set(gca,'ytick',[])
        title('Zoomed in: Change in Orofacial Movements Over Time')
        xlabel('Time (s)')
        
        %add to video file
        drawnow; 
        F(frame) = getframe(gcf); 
        writeVideo(vidfile,F(frame));
        pause(0.1)
    end
    close(vidfile)
end