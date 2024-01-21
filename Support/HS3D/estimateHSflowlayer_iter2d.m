function uv = estimateHSflowlayer_iter2d(frame1,frame2,uv,lambda, maxwarping)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~exist('lambda','var')
    lambda = 80;
end
% lambda = 80;
% lambda = 20;
if ~exist('maxwarping','var')
    maxwarping = 100;
end

H = size(frame1,1);
W = size(frame1,2);

npixels = H*W;
[x, y] = meshgrid(1:W,1:H);

% Kernel to get gradient
h = [1 -8 0 8 -1]/12;

% h = [-1,-2,-1;0,0,0;1,2,1];

% ave_filter = [0,1,0;1,0,1;0,1,0];

a = repmat(lambda,npixels,1);
for i = 1 : maxwarping
    x0 = x+uv(:,:,1);
    y0 = y+uv(:,:,2);
    
    warpimg2 = [];
    warpimg2 = interp2_bicubic(frame2,x0,y0);
  
    warpimg2(isnan(warpimg2)) = 0;
    
    It = warpimg2 - frame1;
    Ix = imfilter(warpimg2,h);
    Iy = imfilter(warpimg2,h.');
    
    Ix = reshape(Ix,npixels,1);
    Iy = reshape(Iy,npixels,1);
    It = reshape(It,npixels,1);
    
    %% average using 5*5 mask
%     u_ave = medfilt2(uv(:,:,1),[5 5]);
%     v_ave = medfilt2(uv(:,:,2),[5 5]);
%     u_ave = imfilter(uv(:,:,1),ave_filter);
%     v_ave = imfilter(uv(:,:,2),ave_filter);

    U = reshape(uv(:,:,1),npixels,1);
    V = reshape(uv(:,:,2),npixels,1);
    
    %%
    U_refine = U - (Ix.*(Ix.*U + Iy.*V + It))./(a.*a + Ix.*Ix + Iy.*Iy);
    V_refine = V - (Iy.*(Ix.*U + Iy.*V + It))./(a.*a + Ix.*Ix + Iy.*Iy);
    
    
    uv(:,:,1) = reshape(U_refine, size(uv(:,:,1)));
    uv(:,:,2) = reshape(V_refine, size(uv(:,:,2)));
    
    uv(:,:,1) = medfilt2(uv(:,:,1),[7 7]);
    uv(:,:,2) = medfilt2(uv(:,:,2),[7 7]);
    
end

end

