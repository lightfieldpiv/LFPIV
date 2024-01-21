% Deform function for NGP
% Input: source location and target location
% output: index_compute
function index_compute = NGP_Deform(s_frame_point,t_frame_point)

%m = size(s_frame_point,1);
% index = randperm(m,m/10);

original = pointCloud(s_frame_point);
% sample   = pointCloud(s_frame_point(index,:));

gridStep = 0.001;

% sample = pcdownsample(original,'random',0.2);
sample = pcdownsample(original,'gridAverage',gridStep);

%figure;pcshow(original);
%figure;pcshow(sample,'MarkerSize',200);

index_tmp_control = knnsearch(s_frame_point, sample.Location,'k',1); % extract index

index_compute = NRG_graph( s_frame_point,t_frame_point,index_tmp_control);

end

