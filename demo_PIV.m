
clc;
clear;
close all;

%% MIAN FUNCTION (PIV)
addpath(genpath('support'));

%% Input point cloud from two frames or gernerate it
folder = './test_data/';
frameNameS = 4447;
frameNameT = 4457;
numFrame = frameNameT - frameNameS +1;
pt_frame = cell(numFrame,1);

for i = 1:numFrame
    pt_frame{i} = pcread([folder num2str(frameNameS+i-1) '_UWN.ply']);
end

%% PIV main function
% the size of the volumn.
sample_ratio = 5;
param.H_size = 80/sample_ratio;
param.W_size = 120/sample_ratio;
param.D_size = 120/sample_ratio;
% resize the volume
frame_pos_resize  = BatchScale(pt_frame, param.H_size, param.W_size, param.D_size);
%% Reconstruct the Dynamic Fluid
% flag = 1;% regular optical flow
% flag = 2;% gt corres optical flow
% flag = 3;% gt corres & div optical flow
% flag = 4;% div optical flow
% flag = 5;% calculate corres optical fdlow
% flag = 6;% calculate corres optical flow & div-free
param.flag = 6;
uv = cell(numFrame-1,1);
for i = 1:numFrame-1
    %   find the point correspondence between two point clouds
    param.PointCorespondenceIndex = NGP_Deform(pt_frame{i}.Location, pt_frame{i+1}.Location);
    
    uv{i} = FluidFlowRecon(frame_pos_resize{i}, frame_pos_resize{i+1}, param);
    % display the result
    plotFlow3(param.W_size, param.H_size, param.D_size, uv{i}(:,:,:,1), uv{i}(:,:,:,2), uv{i}(:,:,:,3), [5,5,5], 1);
    title([num2str(frameNameS) '-' num2str(frameNameT) '-PIV-flag-' num2str(param.flag)]);
end
% visualize the fluid as path line figure
% savePath = './result/';
% SaveFolder = [savePath num2str(frameNameS) '_' num2str(frameNameT) '_PIV_flag' num2str(param.flag) '\'];
% if ~exist(SaveFolder,'dir')
%     mkdir(SaveFolder);
% end
% VisPathLineFluid(uv, param, SaveFolder);


