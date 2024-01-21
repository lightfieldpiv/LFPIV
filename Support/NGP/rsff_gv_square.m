function energy = rsff_gv_square(x,control_point,source_obj,...
                                 target_obj,weight,weightTrans,...
                                 coeff_rigid ,coeff_smooth,alpha_point)
% rsff_gv no linear optimization for non rigid registration
% Evaluate the function.
%tic
% global source_face
% global index_border
%% E_rigid

control_point_num = size(control_point,1);

%% E_rigid
% disp('E_rigid:');
% tic
x_shape = reshape(x,12,size(control_point,1));

a1 = x_shape(1:3,:);
a2 = x_shape(4:6,:);
a3 = x_shape(7:9,:);
b  = x_shape(10:12,:);

one_rigid = ones(1,size(a1,2));

rigid4 = (one_rigid - sum(a1.*a1,1))';
rigid5 = (one_rigid - sum(a2.*a2,1))';
rigid6 = (one_rigid - sum(a3.*a3,1))';

rigid1 = sum(a1.*a2,1)';
rigid2 = sum(a1.*a3,1)';



rigid3 = sum(a2.*a3,1)';

E_rigid_square = sqrt(coeff_rigid)*[rigid1 ;rigid2 ;rigid3 ;rigid4 ;rigid5 ;rigid6];
% t1 = (a1.*a2.^2 ;

%% E_smooth

% disp('r_diance:');
% 
% tic
control_trans = control_point';

E_smooth_square = zeros(control_point_num*control_point_num,3);
smooth_count = 1;
for i = 1 : control_point_num
    for j = 1 : control_point_num
     A_i = [a1(:,i),a2(:,i),a3(:,i)];
     
     x_i  = control_trans(:,i);
     
     x_j = control_trans(:,j);
     
     b_i = b(:,i);
     
     b_j = b(:,j);  
     
     tmp = A_i * (x_j - x_i) + x_i + b_i -(x_j + b_j);
     
  
    tmp_square_x = sqrt(weight(i,j)) * tmp(1);
    tmp_square_y = sqrt(weight(i,j)) * tmp(2);
    tmp_square_z = sqrt(weight(i,j)) * tmp(3);
    
    
    E_smooth_square(smooth_count,:) = sqrt(coeff_smooth) * [tmp_square_x,tmp_square_y,tmp_square_z];
    
    smooth_count = smooth_count + 1;
    end
    % x_j = repmat
   % E_smooth(i,j) = 
end

E_smooth_square =  reshape(E_smooth_square',control_point_num*control_point_num*3,1);

% toc;
%% E_fit term
% disp('r_diance:');
% tic
% global index_c0 index_c1

% alpha_point = 0.1;
% alpha_plane = 1.0;
% E_fit = 0;

v_RT =  ApplyTrans(source_obj,a1,a2,a3,b,control_point,weightTrans);

%% if only use control points?
%v_RT = v_RT(index_control,:);

%%
% v_RT = v_RT(index_c0,:);
% threshold = 0.07;
% compute closed point after apply transformation
[index_closest,D] = knnsearch(target_obj,v_RT,'k',1);

% thres_ind = find(D>threshold);
%% discard some outlier correspondence 
target_closest = target_obj(index_closest,:);
% target_normal_index = target_normal(index_closest,:);
% [source_normal,nu] = compute_normal(v_RT,source_face);
% erase normal dot product sum < 0, indicate wrong direction
% dot_normal = target_normal_index.*source_normal;
% dot_sum = sum(dot_normal,2);
% thres_normal = find(dot_sum<0);
%% discard edge correspondence
% % border_vertex = target_obj(index_border,:);
% [border_erase,D] = knnsearch(target_closest,border_vertex,'k',1);
% border_erase = border_erase(:);
% border_erase = unique(border_erase);



% if(size(thres_normal,1)> 0.7*size(v_RT,1))
%     thres_normal = find(dot_sum>0);
% end
    
diff_v = v_RT - target_closest;
%% erase threshold, border , and normal inconsistent repectively
% diff_v(thres_ind,:) = 0;
% diff_v(border_erase,:) = 0;
% diff_v(thres_normal,:) = 0;
% 
% V_erase = v_RT(border_erase,:);
% T_erase = target_closest(border_erase,:);
% 
% fileID = fopen('V_erase.obj','w');
% for i = 1 : size(V_erase,1)
%     fprintf(fileID,'v %.6f %.6f %.6f %d %d %d\n', V_erase(i,1),V_erase(i,2),V_erase(i,3),...
%                                                 1,0,0);
% end
% fclose(fileID);
% 
% fileID = fopen('T_erase.obj','w');
% for i = 1 : size(T_erase,1)
%     fprintf(fileID,'v %.6f %.6f %.6f %d %d %d\n', T_erase(i,1),T_erase(i,2),T_erase(i,3),...
%                                                 0,1,0);
% end
% fclose(fileID);


%%
E_fit1_square =  sqrt(alpha_point) * reshape(diff_v',size(diff_v,1)*3,1); 


% E_fit2_square =   sqrt(alpha_plane) * abs(sum(target_normal_index.* diff_v,2));


energy = [E_rigid_square;E_smooth_square;E_fit1_square;];

% % E_rigid_square'*E_rigid_square
% % E_smooth_square'*E_smooth_square
% % energy'*energy
% toc;
end





