function [E1] = CostVol(LF,center_view, ...
                        numLabel,...
                        LabelUint,winsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Sc = [center_view,center_view];

[num_s, num_t, row, col, ch] = size(LF);

WindowSize = [winsize,winsize];

Ic = squeeze(LF(center_view,center_view,:,:,1:3));

Ic=min(max(Ic,0),1);

[ss, tt] = ndgrid(1:num_s, 1:num_t);
COMB= [ss(:) tt(:)];
aa = strmatch(Sc,COMB);
COMB(aa,:) = [];

% Number of Labels
% numLabel = 75 ;
numTarget = size(COMB,1) ;
% LabelUint = 0.02; % Sub pixel shift
% WindowSize = [3,3];
tau1 = 0.5;
tau2 = 0.5;
alpha = 0.5;

Rho1 = zeros(row, col, numLabel);

%%
% Gx_c = zeros(row, col, ch); Gy_c = zeros(row, col, ch);
Gx_c = zeros(row, col, ch); Gy_c = zeros(row, col, ch);
for ii = 1 : ch
    [Gx_c(:,:,ii), Gy_c(:,:,ii)] = imgradientxy(Ic(:,:,ii));
end



for n = 1 : numTarget
    Sn = COMB(n,:) ;
    Target = im2double(squeeze(LF(Sn(1),Sn(2),:,:,1:3)));
    Target = min(max(Target,0),1);
    
    [Gx_target(:,:,1), Gy_target(:,:,1)] = imgradientxy(Target(:,:,1));
    [Gx_target(:,:,2), Gy_target(:,:,2)] = imgradientxy(Target(:,:,2));
    [Gx_target(:,:,3), Gy_target(:,:,3)] = imgradientxy(Target(:,:,3));
    
    Ftarget_x = fft2(Gx_target);
    Ftarget_y = fft2(Gy_target);
    
    Ftarget = fft2(Target);
    
    Vn = Sn - Sc;
    
    beta = abs(Vn(1))./(abs(Vn(1))+abs(Vn(2)));
    
     parfor ell = 1 : numLabel
        deltar =  - LabelUint * Vn(1) * (ell);
        deltac =  - LabelUint * Vn(2) * (ell);
        
        delta = [deltar deltac];
        In = fn_SubpixelShift(Ftarget, delta, row, col,1);
        SAD1 = min(sqrt( (In(:,:,1) - Ic(:,:,1)).^2 + (In(:,:,2) -...
              Ic(:,:,2)).^2 + (In(:,:,3) - Ic(:,:,3)).^2 ),tau1);
          
        DSI_SAD1 = fn_PatchSum(SAD1,WindowSize);
         
        % Sum of the Gradient Difference
        
        Gx_n = fn_SubpixelShift(Ftarget_x, delta, row, col,0);
        Gy_n = fn_SubpixelShift(Ftarget_y, delta, row, col,0);
        
        GRAD_X1 = min(sqrt( (Gx_c(:,:,1)-Gx_n(:,:,1)).^2 + (Gx_c(:,:,2)-Gx_n(:,:,2)).^2 + (Gx_c(:,:,3)-Gx_n(:,:,3)).^2),tau2);
        GRAD_Y1 = min(sqrt( (Gy_c(:,:,1)-Gy_n(:,:,1)).^2 + (Gy_c(:,:,2)-Gy_n(:,:,2)).^2 + (Gy_c(:,:,3)-Gy_n(:,:,3)).^2),tau2);
        
        DSI_GRAD_X1 = fn_PatchSum(GRAD_X1,WindowSize);
        DSI_GRAD_Y1 = fn_PatchSum(GRAD_Y1,WindowSize);
        DSI_GRAD1 = beta*DSI_GRAD_X1 + (1-beta)*DSI_GRAD_Y1;
        
          
        Rho1(:,:,ell) = (1-alpha).*DSI_SAD1 + alpha.*DSI_GRAD1;
        
     end
    
    if n == 1
        E1 = Rho1;
    else
        E1 = E1 + Rho1;
    end
    
    disp(sprintf('CostVolume... %d / %d', n,numTarget))
end



end

