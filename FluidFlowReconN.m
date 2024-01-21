function [uv] = FluidFlowReconN(s_frame_point,t_frame_point, param)

ratio =4;
corres_data_folder = '../../NGP/'; 

% H = round(434/ratio);
% W = round(625/ratio);
% D = 35;

H = param.H_size + 1;
W = param.W_size + 1;
D = param.D_size + 1;

s = [ H,W,D];
m = prod(s);

debug_translate = 0;
% debug_translate = 1;

try_corres = 0;
flag = param.flag;
%flag = 1;% regular optical flow
% flag = 2;% gt corres optical flow
% flag = 3;% gt corres & div optical flow
%flag = 4;% div optical flow
% flag = 5;% calculate corres optical flow
% flag = 6;% calculate corres optical flow & div-free
%% downsize data
% s_pro = downsizecloud(s_frame_point,ratio);
% t_pro = downsizecloud(t_frame_point,ratio);
% s_frame_point = s_pro;
% t_frame_point = t_pro;
%% transfer to volume data format

frame_in  = transfer_data(s_frame_point,H,W,D);
frame_out = transfer_data(t_frame_point,H,W,D);

% save('frame_in.mat','frame_in')
% save('frame_out.mat','frame_out')

%% get gt correspondence data
%% normalize to 0 to 255
frame1_3d = zeros(size(frame_in));
frame2_3d = zeros(size(frame_out));


for i = 1 : size(frame_in,3)
        frame1 = frame_in(:,:,i);
        frame2 = frame_out(:,:,i);

        frame1 = mat2gray(frame1);
        frame2 = mat2gray(frame2);
  
%         frame1 = downSample(frame1,ratio);
%         frame2 = downSample(frame2,ratio);
        
        frame1 = frame1*255;
        frame2 = frame2*255;
%         
%         if(i==1)
%             frame1_3d = zeros(size(frame1,1),...
%                               size(frame1,2),...
%                               size(frame_in,3));
%             frame2_3d = zeros(size(frame1,1),...
%                               size(frame1,2),...
%                               size(frame_in,3));
%         end
        
%         figure;imshow(frame1,[]);
%         figure;imshow(frame2,[]);

        frame1_3d(:,:,i) = frame1;
        frame2_3d(:,:,i) = frame2;

        close all;
end
frame_in  = frame1_3d;
frame_out = frame2_3d;
    
% H = size(frame_in,1);    
% W = size(frame_in,2); 
%% resize

%%
if(flag==5 || flag ==6)
    load([corres_data_folder,'index_compute.mat'])
    [C_u,C_v,C_w,p_out] = Corres3d_compute(s_frame_point,t_frame_point,index_compute,H,W,D);
end
%%
% s = [H,W,D];
% m = prod(s);
% numFrame = 2;
% 
% smooth_sigma = 0.8;
% f3 = fspecial3('gaussian', 2*round(1.5*smooth_sigma) +1 );

if(flag == 1) % non
    uv = estimateHSflow3d(frame_in, frame_out);
elseif(flag ==2 || flag==5) % add corres
    uv = estimateHSflow3d_corres(frame_in, frame_out,C_u,C_v,C_w,0);
elseif(flag ==3 || flag==6) % add corres-div-free
    uv = estimateHSflow3d_corres(frame_in, frame_out,C_u,C_v,C_w,1);
elseif(flag ==4) % add div-free
    uv = estimateHSflow3d_div(frame_in, frame_out);
end
% uv(:,:,:,3) = 0;uv
%% error measure
% ratio=8;
% uv = permute( uv, [3 2 1 4 5]);
% uv = permute( uv, [3 2 1 4]);
u = uv(:,:,:,1);
v = uv(:,:,:,2);
w = uv(:,:,:,3);

for i = 1 : size(uv,3)
   u(:,:,i) = u(:,:,i) * size(uv,3) / (size(uv,3)-i+1);
   v(:,:,i) = v(:,:,i) * size(uv,3) / (size(uv,3)-i+1);
   w(:,:,i) = w(:,:,i) * size(uv,3) / (size(uv,3)-i+1);
end

uv_r(:,:,:,1) = u;
uv_r(:,:,:,2) = v;
uv_r(:,:,:,3) = w;

% W = 125; H=86;D=30; 
plotFlow3(W, H, D, uv(:,:,:,1), uv(:,:,:,2), uv(:,:,:,3), [10,10,4], 1);