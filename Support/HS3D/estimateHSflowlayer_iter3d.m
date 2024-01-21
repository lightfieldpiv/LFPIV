function uv = estimateHSflowlayer_iter3d(frame1,frame2,...
                                         uv,lambda,...
                                         maxwarping)
% 3D optical flow using HS

[H W D] = size(frame1);

numFrame = 2;

npixels = H*W*D;

[x, y, z] = meshgrid(1:W,1:H,1:D);

h = [1 -8 0 8 -1]/12;

a = repmat(lambda,npixels,1);


%% 5*5*5 average mask
m_s = 5;
ave_mask = zeros(m_s,m_s,m_s);

for i = 1 : size(ave_mask,3)
    ave_mask(:,:,i) = ones(m_s,m_s)/(m_s^3);
end

%%
for i = 1 : maxwarping
    
    x0 = x + uv(:,:,:,1);
    y0 = y + uv(:,:,:,2);
    z0 = z + uv(:,:,:,3);
    
    warpimg2 = [];
    warpimg2 = interp_valid(frame2,x0,y0,z0,'cubic');
    
    % calculate derivate
    It = reshape( warpimg2 - frame1, npixels, 1 );
    
    Ix = reshape( imfilter(warpimg2,h, 'symmetric'),...
                 npixels, 1 );
             
    Iy = reshape( imfilter(warpimg2,h','symmetric'),npixels,1);
    
    Iz_tmp = imfilter( permute( warpimg2, [3 1 2]), h',  'symmetric');
    Iz_tmp = permute( Iz_tmp, [2 3 1] );
    Iz = reshape(Iz_tmp,npixels,1 );
            
    % reshape parameter
    U = reshape(uv(:,:,:,1),npixels,1);
    V = reshape(uv(:,:,:,2),npixels,1);
    W = reshape(uv(:,:,:,3),npixels,1);
    
    % minimized the euler-lagrange equation
    tmp = (Ix.* U + Iy.*V + Iz.*W + It)./...
          (a.*a + Ix.*Ix + Iy.*Iy + Iz.*Iz);
    
    U_refine = U - Ix.*tmp;
    V_refine = V - Iy.*tmp;
    W_refine = W - Iz.*tmp;
    
    uv(:,:,:,1) = reshape(U_refine, size(uv(:,:,:,1)));
    uv(:,:,:,2) = reshape(V_refine, size(uv(:,:,:,2)));
    uv(:,:,:,3) = reshape(W_refine, size(uv(:,:,:,3)));
    
    uv(:,:,:,1) = imfilter(uv(:,:,:,1),ave_mask);
    uv(:,:,:,2) = imfilter(uv(:,:,:,2),ave_mask);
    uv(:,:,:,3) = imfilter(uv(:,:,:,3),ave_mask);
    
end


end

