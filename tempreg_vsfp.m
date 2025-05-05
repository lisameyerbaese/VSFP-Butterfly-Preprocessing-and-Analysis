function [A_reg, D_reg] = tempreg_vsfp(A_input,D_input,tform_temp)

%% Function for aligning trial data to aligned template image for each mouse
% Uses imwarp and tform previously calculated on base image in
% preProcVSFP10L script

%% Some things to start
[sX,sY,sZ] = size(A_input);
A_reg = zeros(sX,sY,sZ);
D_reg = zeros(sX,sY,sZ);

map_ref = imref2d([sX,sY]);

%% Warp D_input
for i = 1:sZ 
    D_reg(:,:,i) = imwarp(D_input(:,:,i),tform_temp,'nearest','OutputView',map_ref);
end

%% Warp A_input
for i = 1:sZ
    A_reg(:,:,i) = imwarp(A_input(:,:,i),tform_temp,'nearest','OutputView',map_ref);
end

