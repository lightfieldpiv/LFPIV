
%%
% addpath('tools/')
% addpath('../../SIFT3D/SIFT3D 1.4.5/lib/sift3d/wrappers/matlab')
synthetic_data_folder = '../stable fluid/2D example/3Dstablefluid/result/'; 

load([synthetic_data_folder,'s_point.mat'])
load([synthetic_data_folder,'t_point.mat'])
%%
m = size(s_frame_point,1);
index = randperm(m,m/10);

original = pointCloud(s_frame_point);
% sample   = pointCloud(s_frame_point(index,:));

gridStep = 10;

% sample = pcdownsample(original,'random',0.2);
sample = pcdownsample(original,'gridAverage',gridStep);

figure;pcshow(original);
figure;pcshow(sample);

index_tmp_control = knnsearch(s_frame_point, sample.Location,'k',1); % extract index

NRG_graph( s_frame_point,t_frame_point,index_tmp_control);
