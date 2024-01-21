function [C_u,C_v,C_w,p_out] = Corres3d_compute(p,q,corres,H,W,D)

[p_out,uvw_out] = transfer_data_corres_compute(p,q,corres,H,W,D);

[C_u,C_v,C_w] = convertFlow_compute(uvw_out,p_out,H,W,D);

end

function [c_u,c_v,c_w] = convertFlow_compute(uvw,p,H,W,D)
% construct u v vector
    c_u = zeros(H,W,D);
    c_v = zeros(H,W,D);
    c_w = zeros(H,W,D);
    
    for i = 1 : size(p,1)
            c_u(p(i,2),p(i,1),p(i,3)) = uvw(i,1); 
            c_v(p(i,2),p(i,1),p(i,3)) = uvw(i,2);
            c_w(p(i,2),p(i,1),p(i,3)) = uvw(i,3);
    end
    
end

function [p_out,uv_out] = transfer_data_corres_compute(p,q,corres,H,W,D)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here


% frame_out = zeros(H,W,D);

p_out = zeros(size(p,1)*8,3);
uv_out = zeros(size(p,1)*8,3);
%v_out = zeros(size(p,1)*8,3);

count = 1 ;

q_corres = q(corres,:);

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
    corner3d = [j0,i0,k0;... % y , x , z
                j0,i1,k0;...
                j1,i0,k0;...
                j1,i1,k0;...
                j0,i0,k1;...
                j0,i1,k1;...
                j1,i0,k1;...
                j1,i1,k1;];
    
    
    % calculate u and v, attention , u and v first coordinate is column,
    % second is row
%     u_tmp = interp3(u, y,x,z);
%     v_tmp = interp3(v, y,x,z);
    %
    
    u_tmp = q_corres(i,1) - p(i,1); % x
    v_tmp = q_corres(i,2) - p(i,2); % y
    w_tmp = q_corres(i,3) - p(i,3); % z
    
    for j = 1 : 8
              
        p_out(count,:) = [corner3d(j,2),... % x
                          corner3d(j,1),... % y
                          corner3d(j,3)]; % z
                      
%         u_tmp = u(corner3d(j,1),corner3d(j,2),corner3d(j,3));
%         v_tmp = v(corner3d(j,1),corner3d(j,2),corner3d(j,3));
%         w_tmp = w(corner3d(j,1),corner3d(j,2),corner3d(j,3));
        
        uv_out(count,:) = [u_tmp,v_tmp,w_tmp];               
                      
        count = count + 1;
    end
    
%     if(sum(weight_all)> 1.0001)
%           disp('here');
%     end

end
end