function frame_out = transfer_data(p,H,W,D)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% frame_out = zeros(H,W,D);
frame_out = zeros(H,W,D);

for i = 1 : size(p,1)
    x = p(i,1);y=p(i,2);z=p(i,3);
    %% for x
    if(x < 0.5+1)
      x = 0.5+1;
    end
        
    if(x> W+0.5-1)
       x = W+0.5-1;
    end
    %% for y
    if(y < 0.5+1)
       y = 0.5+1;
    end
        
    if(y > H+0.5-1)
       y = H+0.5-1;
    end
%% for Z
   if(z < 0.5+1)
       z = 0.5+1;
   end
        
   if(z > D+0.5-1)
       z = D+0.5-1;
   end
     
    %%
   i0 = floor(x);
   i1 = i0 + 1;
    
    j0 = floor(y);
    j1 = j0 + 1;
    
    k0 = floor(z);
    k1 = k0 + 1;
    
%     corner3d = [i0,j0,k0;...
%                 i1,j0,k0;...
%                 i0,j1,k0;...
%                 i1,j1,k0;...
%                 i0,j0,k1;...
%                 i1,j0,k1;...
%                 i0,j1,k1;...
%                 i1,j1,k1;];
    corner3d = [j0,i0,k0;...
                j0,i1,k0;...
                j1,i0,k0;...
                j1,i1,k0;...
                j0,i0,k1;...
                j0,i1,k1;...
                j1,i0,k1;...
                j1,i1,k1;];
    
    weight = zeros(8,1);
    
    for j = 1 : 8
%         weight(j) = (sqrt(3) - norm(corner3d(j,:)-[x,y,z]))^3;
           weight(j) = (sqrt(3) - norm(corner3d(j,:)-[y,x,z]))^3;
        
        if(weight(j)<0)
            disp('here');
        end
    end
    %
    weight_all = zeros(8,1);
    for j = 1 : 8
        weight_all(j) = weight(j)/sum(weight);
        
        frame_out(corner3d(j,1),...
                  corner3d(j,2),...
                  corner3d(j,3)) = weight_all(j);
    end
    
    if(sum(weight_all)> 1.0001)
          disp('here');
    end

end




end

