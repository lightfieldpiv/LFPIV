function [index, cnt] = ParticleDepth(fileloc,start_slab,end_slab,numLabel,LabelUnit, winsize)

%% hist equal the data
load(fileloc)
LF_hist = LFHistEqualize(LF);
LF_center = LF_hist(start_slab:end_slab,start_slab:end_slab,:,:,1:3);
% LF_center = DownSampleLF(LF,ratio,start_slab,end_slab);
center_view = (end_slab - start_slab)/2 + 1;

%% detect particle location
img_center = squeeze(LF_center(center_view,center_view,:,:,:));
% %% shift the image if necessary
% LF_shifted = ShiftLF(LF_center,...
%                       start_slab,end_slab,...
%                       shift_pixel);
% LF_center = LF_shifted;

%% detect particle
thres = 30;
cnt = detectParticleLoc(img_center,thres);

%% estimate depth cost volume
% winsize = 1;
tic;
[E1] = CostVol(LF_center,center_view,numLabel,LabelUnit,winsize);
toc;
%% aggregation
% Ic = im2double(squeeze(LF_center(center_view,center_view,:,:,1:3))); % Guided image 
% param.r = end_slab - start_slab+1;% illum for 7, ligtro for 5
% param.eps = 0.0001;
% tic;
% E2 = CostAgg(E1,Ic,param);
% toc;
%% find the depth label with the lowest cost
[test, index] = findLoweset(E1);
% [test, index] = findLoweset(E2);
figure;imshow(index,[]);
radi = 3*ones(size(cnt,1),1);
h = viscircles(cnt(:,1:2),radi);
figure;imshow(img_center,[]);
h = viscircles(cnt(:,1:2),radi);
figure;imshow(img_center,[]);
figure;imshow(index,[]);



