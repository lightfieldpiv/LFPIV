function uv = estimateHSflow3d(frame_source,frame_target)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
[H,W,D] = size(frame_source);

pyramid_level = 3;
smooth_sigma = 0.8;

numFrame = 2;

f3 = fspecial3('gaussian', 2*round(1.5*smooth_sigma)+1);

pyramid1 = cell(pyramid_level,1);
pyramid2 = cell(pyramid_level,1);

pyramid1{1} = frame_source;
pyramid2{1} = frame_target;

pyramid_size = zeros(pyramid_level,3);
pyramid_size(1,:) = size(frame_source);

for m = 2: pyramid_level
    pyramid_size(m,:) = pyramid_size(m-1,:).*[0.5 0.5 0.5];
    
    pyramid1{m-1} = imfilter(pyramid1{m-1},f3);
    pyramid2{m-1} = imfilter(pyramid2{m-1},f3);
    
    pyramid1{m} = resize(pyramid1{m-1},pyramid_size(m,:));
    pyramid2{m} = resize(pyramid2{m-1},pyramid_size(m,:));
    
    %figure;imshow(pyramid1{m-1}(:,:,1),[])
    %figure; imshow(pyramid2{m-1}(:,:,1),[])
end

% coarst-to-fine compute the flow
uv = zeros(H,W,D,3);
lambda = 100;
for levels = pyramid_level:-1:1
    fprintf('level: %d \n',levels);
    uv = resample_flow3d( uv, pyramid_size(levels,:));
    
    uv = estimateHSflowlayer3d(pyramid1{levels},pyramid2{levels},uv,lambda);
end

end

