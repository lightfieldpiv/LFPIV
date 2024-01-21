function uv = estimateHSflowlayer3d(source,target,uv,lambda)
% maxwarping

[H,w,D] = size(source);

npixels = H*w*D;

[x,y,z] = meshgrid(1:w,1:H,1:D);

% build differential matrix and laplacian matrix according to image size
e = ones(npixels,1);

dy = spdiags([-e e],0:1,npixels,npixels);
dx = spdiags([-e e],[0, H],npixels,npixels);
dz = spdiags([-e e],[0, H*w],npixels,npixels);

dy(H:H:npixels,:) = 0;
for i=1:D
   dx( H*(w*i-1)+1 : 1 : H*w*i, : ) = 0; 
end
dz( H*w*(D-1)+1:1:H*w*D , : ) = 0;
% L = dx.'*dx + dy.'*dy + dz.'*dz;
L = dx'*dx + dy'*dy + dz'*dz;

% Kernel to get gradient
h = [1 -8 0 8 -1]/12;
% h = [1 -1];

%% 5*5*5 average mask
m_s = 5;
ave_mask = zeros(m_s,m_s,m_s);

for i = 1 : size(ave_mask,3)
    ave_mask(:,:,i) = ones(m_s,m_s)/(m_s^3);
end

ave_mask = ones(m_s,m_s,3)/(m_s*m_s*3);

for i=1:3
     x1 = x + uv(:,:,:,1);
     y1 = y + uv(:,:,:,2);
     z1 = z + uv(:,:,:,3);
     warpimg2 = [];
     warpimg2 = interp_valid(target,x1,y1,z1,'cubic');
     
      % TODO#4: compute image gradient Ix, Iy, and Iz 
    It = reshape( warpimg2 - source, npixels, 1 );
    Ix = spdiags( reshape( imfilter(warpimg2,h, 'symmetric'), npixels, 1 ),0,npixels,npixels);
    Iy = spdiags( reshape( imfilter(warpimg2,h','symmetric'), npixels, 1 ),0,npixels,npixels);
    Iz_tmp = imfilter( permute( warpimg2, [3 1 2]), h',  'symmetric');
    Iz_tmp = permute( Iz_tmp, [2 3 1] );
    Iz = spdiags( reshape(Iz_tmp,npixels,1 ),0,npixels,npixels);
    
     U = reshape(uv(:,:,:,1),npixels,1);
     V = reshape(uv(:,:,:,2),npixels,1);
     W = reshape(uv(:,:,:,3),npixels,1);
    
    A = [Ix*Ix+lambda*L Ix*Iy Ix*Iz;...
        Ix*Iy Iy*Iy+lambda*L Iy*Iz;...
        Ix*Iz Iy*Iz Iz*Iz+lambda*L];
    
    b = -[Ix*It+lambda*L*U;...
          Iy*It+lambda*L*V;...
          Iz*It+lambda*L*W;];
    
    L1 = ichol(A);
    [deltauv,~] = pcg(A, b, 1e-4, 100, L1, L1');  
%     deltauv = A\b;
    
    deltauv = reshape( deltauv, size(uv));
    
    deltauv(deltauv>1) = 1;
    deltauv(deltauv<-1) = -1;
    
    uv(:,:,:,1) = uv(:,:,:,1) + deltauv(:,:,:,1);
    uv(:,:,:,2) = uv(:,:,:,2) + deltauv(:,:,:,2);
    uv(:,:,:,3) = uv(:,:,:,3) + deltauv(:,:,:,3);
    
    
     % TODO#6: use median filter to smooth the flow map 
    uv(:,:,:,1) = imfilter(uv(:,:,:,1),ave_mask,'symmetric');
    uv(:,:,:,2) = imfilter(uv(:,:,:,2),ave_mask,'symmetric');
    uv(:,:,:,3) = imfilter(uv(:,:,:,3),ave_mask,'symmetric');
%     B = medfilt3(A,[m n p])
    
%     uv(:,:,:,1) = medfilt3(uv(:,:,:,1),[5 5 3]);
%     uv(:,:,:,2) = medfilt3(uv(:,:,:,2),[5 5 3]);
%     uv(:,:,:,3) = medfilt3(uv(:,:,:,3),[5 5 3]);
%     
    fprintf('Warping step: %d, Incremental u norm: %3.5f \n', i,   sum(sum(sum(abs(deltauv(:,:,:,1))))));
    fprintf('Warping step: %d, Incremental v norm: %3.5f \n', i,   sum(sum(sum(abs(deltauv(:,:,:,2))))));
    fprintf('Warping step: %d, Incremental w norm: %3.5f \n\n', i, sum(sum(sum(abs(deltauv(:,:,:,3))))));
end

end