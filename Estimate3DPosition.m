function [S_location,S_disparity, cnt] = Estimate3DPosition(Folder, img_No, K)
%% This Script is used to estimate the 3D position of particles(Zhong Li's paper)

% Folder1 = 'D:/Lytro2Data/Images3/AAA304/Test_TB_UW_P/';
img_name = ['IMG_' img_No '__Decoded.mat'];

%% some superparameters
% the index of subaperture images which is used
start_slab = 5;
end_slab = 11;
winsize = 2;
LabelUnit = 0.02;
% depth resolution
numLabel = 50;

%% compute the disparity Particle
[S_disparity, cnt] = ParticleDepth([Folder, img_name],start_slab,end_slab,numLabel,LabelUnit,winsize);

%% compute the 3D point by the intrinsic parameters
S_location =  DisToDepthCAM(K, S_disparity, cnt, LabelUnit);                       

% %% scale the point cloud with the volume size
% % min x 
% min_x1 = min(S_location1(:,1));
% min_x2 = min(S_location2(:,1));
% min_x  = min(min_x1, min_x2); 
% 
% max_x1 = max(S_location1(:,1));
% max_x2 = max(S_location2(:,1));
% max_x  = max(max_x1, max_x2);
% 
% % min y
% min_y1 = min(S_location1(:,2));
% min_y2 = min(S_location2(:,2));
% min_y  = min(min_y1, min_y2); 
% 
% max_y1 = max(S_location1(:,2));
% max_y2 = max(S_location2(:,2));
% max_y  = max(max_y1, max_y2);
% 
% % min z
% min_z1 = min(S_location1(:,3));
% min_z2 = min(S_location2(:,3));
% min_z  = min(min_z1, min_z2); 
% 
% max_z1 = max(S_location1(:,3));
% max_z2 = max(S_location2(:,3));
% max_z  = max(max_z1, max_z2);

% the volume size
% H_size = 100;
% W_size = 125;
% D_size = 20;
% s1_resize = scalePoints(S_location1,min_x,max_x,min_y,max_y,min_z,max_z,H_size,W_size,D_size);                        
% s2_resize = scalePoints(S_location2,min_x,max_x,min_y,max_y,min_z,max_z,H_size,W_size,D_size);           
% 
% % display the resize 3D particle.
% figure;scatter3(s1_resize(:,1), s1_resize(:,2), s1_resize(:,3), 10, s1_resize(:,3), 'filled');
% figure;scatter3(s2_resize(:,1), s2_resize(:,2), s2_resize(:,3), 10, s2_resize(:,3), 'filled');
% 
% %ratio = 10;
% % figure;scatter3(S_location1(:,1)/ratio, S_location1(:,2)/ratio, S_location1(:,3), 100, S_location1(:,3), 'filled');% 10 is particle size
% % figure;scatter3(S_location2(:,1)/ratio, S_location2(:,2)/ratio, S_location2(:,3), 100, S_location2(:,3), 'filled');% 