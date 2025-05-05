% this function loads in the preProcessed imaging data and adds to the same
% .mat files the following things 
%   1) PC1 and PC2 of the raw imaging data
%   2) global signal
%   3) low pass filtered hemo signal
%   4) behvaioral data that was collected (if applicable)
%
%   Function Inputs: 'resting', 'ephys' or 'airPuff'

function mergeFile(type)
    
    disp(['** Processing ' type ' data**'])
    %type = 'ephys';
    %type = 'resting';
    %type = 'airPuff';
    [~, name] = system('hostname');
    if contains(name,'jaeger')
        startFile = 'X:\labs\keilholz-lab\Lisa';s
    else
        startFile = 'X:\keilholz-lab\Lisa';
    end
    

    % select image used for mask
    disp('Selecting Cortical Map to Use for Masking')
    path_map = 'C:\Users\lmeyerb\OneDrive - Emory University\Desktop\Allen Atlas';
    file = 'InkedmaskCortical_Map.jpg';
    [maskImg,~] = imread(fullfile(path_map, file));
    
    % load in cortical map Cortical_Map.png
    disp('Selecting Cortical Map to Use for Registration')
    file = 'Cortical_Map.png';
    [fixed_map,~] = imread(fullfile(path_map, file));

    if isequal(type, 'resting')
        addpath([startFile, '\VSFP ButterFly\Info'])
        T = readtable('VSFP_50Hz_Proc.xlsx');
        all_mice = unique(T.Mouse); 
        all_mice(2) = []; % for now delete VSFP 25, no pupil data

        for m = 1: length(all_mice) 
            mouse = all_mice{m};
            ind = find(contains(T.Mouse,mouse));
            all_trials = T.Trials(ind);
            all_dates = T.Date(ind);
        
            for k = 1:length(all_trials)
                %get corresponding pupil data
                FindTable = T((T.Trials == all_trials(k) & (T.Date == all_dates(k))),:);
                date = FindTable.Date;
                trial = FindTable.Trials;
                pupil_file = FindTable.video_proc_file{1,1};
    
                % load pupil data
                saved_data = strcat([startFile, '\VSFP ButterFly\Data\Pupil Diameter\',pupil_file]);
                pupil_data = load(saved_data).outData;
                % get pupil data
                frame_inds= FindTable.start_frame:FindTable.end_frame;
                pupil_data_crop = medfilt1(pupil_data.pupil(frame_inds),10); %apply 10th order filter
                %Upsample pupil data x2
                pupil_data2x = interp(pupil_data_crop,2);
    
                % load face state data
                saved_data = strcat([startFile, '\VSFP ButterFly\Data\Face States\',pupil_file]);
                face_data = load(saved_data).faceData;
                % get pupil data
                frame_inds= FindTable.start_frame:FindTable.end_frame;
                if max(frame_inds) > length(face_data.mvmt_filt_z)
                    face_data_crop = face_data.mvmt_filt_z(min(frame_inds):end);
                else
                    face_data_crop = face_data.mvmt_filt_z(frame_inds);
                end
                face_data2x = interp(face_data_crop,2);
    
                % load the raw imaging data
                image_file = strcat(mouse,'_',num2str(date),'_',num2str(trial));
                image_data = strcat(startFile, '\VSFP ButterFly\Data\VSFP_50Hz\',image_file,'.mat');
                load(image_data);
                rawDonor = out.imgD;
                rawAcceptor = out.imgA;
                fs = out.Fs;

                % get the voltage and hemo singal using the projection
                % method
                [projectedHemo, projectedVolt] = projectionMethod(rawDonor.* out.mask, rawAcceptor.* out.mask, fs, image_file, startFile);
                
                % calculate and compare the global signal between ratio metric and
                % projection method
                gsHemoProjected = calcGlobalSingal(projectedHemo);
                gsVoltProjected = calcGlobalSingal(projectedVolt);
                
                voltData = out.imgDR3;
                hemoData = out.imgDS;
                gsVolt = calcGlobalSingal(voltData .* out.mask);
                gsHemo = calcGlobalSingal(hemoData.* out.mask);

                f1 = figure(1);
                subplot(2,1,1)
                plot(zscore(gsVoltProjected), 'DisplayName', 'Volt Projected')
                hold on
                plot(zscore(gsVolt),'DisplayName', 'Volt Ratiometric')
                title('Voltage')
                legend()
                subplot(2,1,2)
                plot(zscore(gsHemoProjected), 'DisplayName', 'Hemo Projected')
                hold on
                plot(zscore(gsHemo),'DisplayName', 'Hemo Ratiometric')
                legend()
                title('Hemo')
                close(f1)
                    
                % low pass filter the hemo signal to 5Hz 
                hemo2D = reshape(hemoData, 10000,[]); 
                lpHemo2D = lowpass(hemo2D',5,fs); %filters each column independantly
                hemoLP = reshape(lpHemo2D',100,100,[]);
                gsHemoLP = lowpass(gsHemo',5,fs);

                % low pass filter the hemo signal to 1Hz 
                hemo2D = reshape(hemoData, 10000,[]); 
                lpHemo2D = lowpass(hemo2D',1,fs); %filters each column independantly
                hemoLP_1Hz = reshape(lpHemo2D',100,100,[]);
                gsHemoLP_1Hz = lowpass(gsHemo',1,fs);

                % add all the variables to the saved .mat file
                out.hemoLP = hemoLP;
                out.gsVolt = gsVolt;
                out.gsHemo = gsHemo;
                out.gsHemoLP = gsHemoLP;
                out.projectedHemo = projectedHemo;
                out.projectedVolt = projectedVolt;
                out.hemoLP1Hz = hemoLP_1Hz;
                out.gsHemoLP1Hz = gsHemoLP_1Hz;

                % save all the new variables to the same data file
                save(image_data,'out','pupil_data2x', 'face_data2x','-v7.3')

            end
        end

    elseif isequal(type,'airPuff')

        allTrials = dir([startFile,'\VSFP ButterFly\Data\VSFP_WhiskerStim\VSFP00_0917\']);
        %load in information about the stimulus type and stimulus side
        for i = 3:length(allTrials)
            load([allTrials(i).folder,'\',allTrials(i).name])
            
            rawDonor = out.imgD;
            rawAcceptor = out.imgA;
            fs = out.Fs;
            [projectedHemo, projectedVolt] = projectionMethod(rawDonor.* out.mask, rawAcceptor.* out.mask, fs, allTrials(i).name(1:end-4), startFile);
            
            % compare the global signal between ratio metric and
            % projection method
            gsHemoProjected = calcGlobalSingal(projectedHemo);
            gsVoltProjected = calcGlobalSingal(projectedVolt);

            voltData = out.imgDR3;
            hemoData = out.imgDS;
            gsVolt = calcGlobalSingal( voltData .* out.mask);
            gsHemo = calcGlobalSingal( hemoData.* out.mask);

            figure()
            subplot(2,1,1)
            plot(zscore(gsVoltProjected), 'DisplayName', 'Volt Projected')
            hold on
            plot(zscore(gsVolt),'DisplayName', 'Volt Ratiometric')
            title('Voltage')
            legend()
            subplot(2,1,2)
            plot(zscore(gsHemoProjected), 'DisplayName', 'Hemo Projected')
            hold on
            plot(zscore(gsHemo),'DisplayName', 'Hemo Ratiometric')
            title('Hemo')
            legend()
            
            % low pass filter the hemo signal to 5Hz 
            hemo2D = reshape(hemoData, 10000,[]); 
            lpHemo2D = lowpass(hemo2D',5,fs); %filters each column independantly
            hemoLP = reshape(lpHemo2D',100,100,[]);
            gsHemoLP = lowpass(gsHemo',5,fs);

            % low pass filter the hemo signal to 1Hz 
            hemo2D = reshape(hemoData, 10000,[]); 
            lpHemo2D = lowpass(hemo2D',1,fs); %filters each column independantly
            hemoLP_1Hz = reshape(lpHemo2D',100,100,[]);
            gsHemoLP_1Hz = lowpass(gsHemo',1,fs);

            % add all the variables to the saved .mat file
            out.hemoLP = hemoLP;
            out.gsVolt = gsVolt;
            out.gsHemo = gsHemo;
            out.gsHemoLP = gsHemoLP;
            out.projectedHemo = projectedHemo;
            out.projectedVolt = projectedVolt;
            out.hemoLP1Hz = hemoLP_1Hz;
            out.gsHemoLP1Hz = gsHemoLP_1Hz;

            % save all the new variables to the same data file
            save([allTrials(i).folder,'\',allTrials(i).name],'out','-v7.3')
            close all
        end

    elseif isequal(type,'ephys')
        allTrials = dir([startFile,'\VSFP ButterFly\Data\VSFP_Ephys\']);
        for i = 3:length(allTrials)
            load([allTrials(i).folder,'\',allTrials(i).name])
            
            rawDonor = out.imgD;
            rawAcceptor = out.imgA;
            fs = out.Fs;
            [projectedHemo, projectedVolt] = projectionMethod(rawDonor.* out.mask, rawAcceptor.* out.mask, fs, allTrials(i).name(1:end-4), startFile);
            
            % compare the global signal between ratio metric and
            % projection method
            gsHemoProjected = calcGlobalSingal(projectedHemo);
            gsVoltProjected = calcGlobalSingal(projectedVolt);

            voltData = out.imgDR3;
            hemoData = out.imgDS;
            gsVolt = calcGlobalSingal( voltData .* out.mask);
            gsHemo = calcGlobalSingal( hemoData.* out.mask);

            % low pass filter the hemo signal to 5Hz 
            hemo2D = reshape(hemoData, 10000,[]); 
            lpHemo2D = lowpass(hemo2D',5,fs); %filters each column independantly
            hemoLP = reshape(lpHemo2D',100,100,[]);
            gsHemoLP = lowpass(gsHemo',5,fs);

            % low pass filter the hemo signal to 1Hz 
            hemo2D = reshape(hemoData, 10000,[]); 
            lpHemo2D = lowpass(hemo2D',1,fs); %filters each column independantly
            hemoLP_1Hz = reshape(lpHemo2D',100,100,[]);
            gsHemoLP_1Hz = lowpass(gsHemo',1,fs);

            % add all the variables to the saved .mat file
            out.hemoLP = hemoLP;
            out.gsVolt = gsVolt;
            out.gsHemo = gsHemo;
            out.gsHemoLP = gsHemoLP;
            out.projectedHemo = projectedHemo;
            out.projectedVolt = projectedVolt;
            out.hemoLP1Hz = hemoLP_1Hz;
            out.gsHemoLP1Hz = gsHemoLP_1Hz;

            % save all the new variables to the same data file
            save([allTrials(i).folder,'\',allTrials(i).name],'out','sortedLFP','-v7.3')
        end

    end

end

function [projectedHemo, projectedVolt] = projectionMethod(rawDonor, rawAcceptor, fs, fileName, startFile)
    % this function takes data from the raw donor and acceptor channels,
    % and then calcualtes the voltage and hemodynamic signal using the
    % projection method outlined in [Carandini et al., 2015]
    
    % computed the relative intensity of the acceptor and donor channels,
    % using first 200ms as the baseline

    [xDim, yDim, zDim] = size(rawDonor);
    avgingFr = floor(fs./5);
    donor2D = reshape(rawDonor, xDim * yDim, zDim);
    baselineDonor = mean(donor2D(:,1:avgingFr),2);
    dfDonor = (donor2D - baselineDonor) ./ baselineDonor;
    dfDonor1D = reshape(dfDonor, 1, []);

    acceptor2D = reshape(rawAcceptor, xDim * yDim, zDim);
    baselineAcceptor = mean(acceptor2D(:,1:avgingFr),2);
    dfAcceptor = (acceptor2D - baselineAcceptor)./ baselineAcceptor;
    dfAcceptor1D = reshape(dfAcceptor, 1, []);

    combinedMatrix = [dfDonor1D ; dfAcceptor1D];

    f1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,3,1)
    s = scatter(dfAcceptor1D,dfDonor1D,'filled','k');
    s.AlphaData = 0.5;
    axis equal
    axis square
    xlabel('Acceptor dF/F (%)')
    ylabel('Donor dF/F (%)')
    title('Scatter Plot of Raw Data')

    % run PCA on this 2D matrix 
    [coeffPCA,scorePCA,~,~,~] = pca(combinedMatrix');
    
    % visualize the projected data
    subplot(1,3,2)
    s = scatter(scorePCA(:,1),scorePCA(:,2),'filled','k');
    s.AlphaData = 0.5;
    axis equal
    axis square
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    title('PCs Volt')

    % visualize the new recasted data along the PC axes 
    reCastData = coeffPCA' * combinedMatrix;
    subplot(1,3,3)
    s = scatter(reCastData(1,:),reCastData(2,:),'filled','k');
    s.AlphaData = 0.5;
    axis equal
    axis square
    title('Recasted Data')

    % reshape the data to make it 3D again
    tempHemo = reCastData(1,:);
    projectedHemo =  reshape(tempHemo, xDim, yDim, zDim);
    tempVolt = reCastData(2,:);
    projectedVolt =  reshape(tempVolt, xDim, yDim, zDim);

    % save all the figures and then close them
    FolderName = [startFile, '\VSFP ButterFly\Data\VSFP_ProjectionMethod\'];
    saveFig(f1, fileName, FolderName)
    close(f1)

end

function [globalSignal] = calcGlobalSingal(data3D)
    [xDim, yDim, zDim] = size(data3D);
    data3D(isnan(data3D))=0;
    data2D = reshape(data3D, xDim * yDim, zDim);
    data2D = zscore(data2D, 0, 2);
    globalSignal = nanmean(data2D, 1);
end