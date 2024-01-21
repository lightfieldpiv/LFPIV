function frame_pos_resize  = BatchScale(frames,H_size,W_size,D_size)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
num = numel(frames);
min_x = zeros(num,1);
min_y = zeros(num,1);
min_z = zeros(num,1);
max_x = zeros(num,1);
max_y = zeros(num,1);
max_z = zeros(num,1);
for  i = 1:numel(frames)
    min_x(i) = frames{i}.XLimits(1);
    max_x(i) = frames{i}.XLimits(2);
    
    min_y(i) = frames{i}.YLimits(1);
    max_y(i) = frames{i}.YLimits(2);
    
    min_z(i) = frames{i}.ZLimits(1);
    max_z(i) = frames{i}.ZLimits(2);
end


%% scale the data
frame_pos_resize = cell(1,num);
for  i = 1:numel(frames)
    frame_pos_resize{i} = scalePoints(frames{i}.Location,...
        min(min_x),max(max_x),...
        min(min_y),max(max_y),...
        min(min_z),max(max_z),...
        H_size,...
        W_size,...
        D_size);
end