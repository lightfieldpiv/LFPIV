function Points = DisToDepthCAM(K, E2, cnt, delta)
% the number of the particle
num = size(cnt,1);
Points = zeros(num,3);
count = 1;
center = [8,8];
for i = 1 : num
        v = floor(cnt(i,2));
        u = floor(cnt(i,1));
        d = E2(v,u)*delta;
        [s1, t1, d1, x1, y1] = RefractRay(K,center(1),center(1),u,v);
        [s2, t2, d2, x2, y2] = RefractRay(K,center(1)+1,center(1),u-d, v);
        
        vp1 = [s1;t1;d1]; dir1 = [x1;y1;1];
        vp2 = [s2;t2;d2]; dir2 = [x2;y2;1];
        A = [dir1,-dir2];
        B = vp2 - vp1;
        ll = pinv(A)*B;
        p1 = vp1 + ll(1)*dir1;
        p2 = vp2 + ll(2)*dir2;
        dist = sqrt(sum((p1-p2).^2));
        
        if abs(dist)<0.001
            X = (p1(1) + p2(1))/2; 
            Y = (p1(2) + p2(2))/2;
            Z = (p1(3) + p2(3))/2;
            Points(count,:) = [X;Y;Z];
        else
            Points(count,:) = [nan;nan;nan];
        end
        count = count+1;
end


%% scale the point cloud to a fixed volume H = 100; W = 120; D =20;
% x side
% s_resize = zeros(size(s_pos));
% 
% H_size = 100;
% W_size = 120;
% D_size = 20;
% 
% min_x = min(s_pos(:,1));
% max_x = max(s_pos(:,1));
% 
% min_y = min(s_pos(:,2));
% max_y = max(s_pos(:,2));
% 
% min_z = min(s_pos(:,3));
% max_z = max(s_pos(:,3));
% 
% for i = 1 : size(s_pos,1)
%     s_resize(i,1) = (s_pos(i,1) - min_x)/(max_x - min_x) * W_size + 1;
%     s_resize(i,2) = (s_pos(i,2) - min_y)/(max_y - min_y) * H_size + 1;
%     s_resize(i,3) = (s_pos(i,3) - min_z)/(max_z - min_z) * D_size + 1;
% 
% end

end


