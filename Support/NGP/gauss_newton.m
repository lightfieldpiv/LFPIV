%function energy = rsff_gv(x)
function [x_out] = gauss_newton(x,control_point,source_obj,target_obj,...
                               weight,weightTrans,coeff_rigid ,coeff_smooth,alpha_point)


%  func = @rsff_gv_square;

initial_sqrt = rsff_gv_square(x,control_point,source_obj,target_obj,weight,weightTrans,coeff_rigid ,coeff_smooth,alpha_point);
  
energy_previous = initial_sqrt'*initial_sqrt;
ini_text = sprintf('initial cost is %f\n',energy_previous);
disp(ini_text);
 
x_shape = reshape(x,12,size(control_point,1));
a1 = x_shape(1:3,:);
a2 = x_shape(4:6,:);
a3 = x_shape(7:9,:);
b  = x_shape(10:12,:);

%--- set fix relection 
%%
for iter = 1 : 20
       
 disp('estimate jacobia matrix....'); 
%   [jac,err] = jacobianest(@(xk)func(xk,control_point,source_obj,...
%              target_obj,target_normal,weight,weightTrans,coeff_rigid ,coeff_smooth,Mdl),x);  
%----------------------------------
v_RT =  ApplyTrans(source_obj,a1,a2,a3,b,control_point,weightTrans);

[index_closest] = knnsearch(target_obj,v_RT,'k',1);

target_closest = target_obj(index_closest,:);

diff_v = v_RT - target_closest;

% symbol = sum(target_normal_index.* diff_v,2);


tic; 
%----------------------------------
jac = Jacobia_C_ready( x,control_point,source_obj,...
                          weight,weightTrans,...
                coeff_rigid ,coeff_smooth,alpha_point);

toc;  
disp('estimate jacobia matrix done!');
   
 
%  jac = sparse(jac_t);
 
%---
tic;  
 disp('solve least square.....');
 F = rsff_gv_square(x,control_point,source_obj,target_obj,weight,weightTrans,coeff_rigid ,coeff_smooth,alpha_point);
%  F_debug = importdata('D:/Zhong_2017_new/NRR/NRR/F.txt');
 
 A = jac'*jac;
 b_vec = -jac'*F;
 %
 %tic;  
 x_delta = A\b_vec;
 disp('solve complete!');
 toc
 
 x_output = x+x_delta';

 F_out = rsff_gv_square(x_output,control_point,source_obj,target_obj,weight,weightTrans,coeff_rigid ,coeff_smooth,alpha_point );
 energy_out = F_out'*F_out;
 ini_text = sprintf('coef_rigid is %f, after num %d iteration cost is %f\n',coeff_rigid,iter,energy_out);
 
 %-------------- adjust ming --------
 compare_value = abs(energy_out - energy_previous)/energy_previous;
 energy_previous = energy_out;
 
 ini_text2 = sprintf('coef_rigid is %f, after num %d iteration diff is %f%\n',coeff_rigid,iter,compare_value);
 
 disp(ini_text);
 disp(ini_text2);
 
  if compare_value < 0.005 && compare_value > 0.000000000001
    disp('stop because threshold set\n');
    break;
 end
 
 x = x_output;
 
%%

%% test result after each iteration
x_shape = reshape(x_output,12,size(control_point,1));

a1 = x_shape(1:3,:);
a2 = x_shape(4:6,:);
a3 = x_shape(7:9,:);
b  = x_shape(10:12,:);
v_RT =  ApplyTrans(source_obj,a1,a2,a3,b,control_point,weightTrans);
%% output result and intermediate result real
%
mkdir('intermedia_result/non_rigid_iteration/');
% transformed point
result_name_source = strcat('intermedia_result/non_rigid_iteration/result_s',num2str(iter-1),'.obj');
result_name_target = strcat('intermedia_result/non_rigid_iteration/result_t',num2str(iter-1),'.obj');
%   options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

fileID = fopen(result_name_source,'w');

for i = 1 : size(v_RT,1)
    fprintf(fileID,'v %.6f %.6f %.6f 255 0 0\n', v_RT(i,1),  v_RT(i,2), v_RT(i,3));
end

fclose(fileID);

fileID = fopen(result_name_target,'w');
for i = 1 : size(target_obj,1)
     fprintf(fileID,'v %.6f %.6f %.6f 128 128 128\n', target_obj(i,1),  target_obj(i,2), target_obj(i,3) );
end


fclose(fileID);
end
x_out = x_output;
 
% t_ori = abs(source_obj-target_obj);
% o_sum = sum(t_ori');
% 
% t_fin = abs(v_RT-target_obj);
% f_sum = sum(t_fin');


end




