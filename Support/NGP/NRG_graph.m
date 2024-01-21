function index_compute = NRG_graph( source,target,index_tmp_control)
%% down_sample

%%
control_point = source(index_tmp_control,:);
%%
index_control = knnsearch(source, control_point,'k',1);

% print source point cloud points and control points number
verbose_source   = sprintf('source points has %d\n',size(source,1));
verbose_control  = sprintf('control points has %d\n',size(control_point,1));

disp(verbose_source);
disp(verbose_control);
%  
%% calculate Radius and weight function for E_smooth term
% global weight
nearest_num = 5;

[~,r_distance] = knnsearch(control_point,source,'k',nearest_num);
[~,r_distance_node] = knnsearch(control_point,control_point,'k',nearest_num);
r_distance = r_distance(:,nearest_num);
r_distance_node = r_distance_node(:,nearest_num);
%-------------------------------------------------------------%
weight = Weight_smooth_ajacent( control_point,r_distance_node);


weightTrans = WeightFunc(source, control_point, r_distance,index_control); % overrall transfomr


%% Compute coefficient for E_fit term
% global index_closest
%% global normal
%% non
num_unkonwn = 12;

x = zeros(num_unkonwn*floor(size(control_point,1)),1);

for i = 1 : length(x)
    if(mod(i,12) == 1)
          x(i) = 1;
     elseif(mod(i,12) == 3)
           x(i) = 0;
           
     elseif(mod(i,12) == 4)
          x(i) = 0;         
     elseif(mod(i,12) == 5)
          x(i) = 1;
          
      elseif(mod(i,12) == 8)
           x(i) = 0;
     elseif(mod(i,12) == 9)
           x(i) = 1;
    elseif(mod(i,12) == 10)
          x(i) = 0;
     elseif(mod(i,12) == 11)
            x(i) = 0;
     elseif(mod(i,12) == 12)
         x(i) = 0;
    else
        x(i)=0;
     end
% x(i)=1;
end

x = x';
%% optimization
% global coeff_rigid coeff_smooth history
coeff_rigid_tmp = 100;
coeff_smooth_tmp = 50;
alpha_point = 1.0;
iter = 1;

initial_sqrt = rsff_gv_square(x,control_point,source,target,weight,weightTrans,coeff_rigid_tmp ,coeff_smooth_tmp,alpha_point);
 
while(1)
    
    coeff_rigid  = coeff_rigid_tmp *(0.5)^(iter-1);
    coeff_smooth = coeff_smooth_tmp *(0.5)^(iter-1);
%     alpha_point = 1.0;
    if coeff_rigid <=49
        break;
    end
   
    x_out =  gauss_newton(x,control_point,source,target,weight,weightTrans,...
                          coeff_rigid,coeff_smooth,alpha_point);

    x = x_out;
     
    iter = iter + 1;
end

x_shape = reshape(x_out,12,size(control_point,1));

a1 = x_shape(1:3,:);
a2 = x_shape(4:6,:);
a3 = x_shape(7:9,:);
b  = x_shape(10:12,:);
v_RT =  ApplyTrans(source,a1,a2,a3,b,control_point,weightTrans);

%% computed

index_compute = knnsearch(target,v_RT,'k',1);

% index_gt = 1:size(target,1);

% comp_error = index_compute - index_gt';

% index_count = find(comp_error==0);
% 
% rate = length(index_count);

% fprintf('accurate rate is %d percent\n',100*rate/size(index_compute,1));
save('index_compute.mat','index_compute');
%% orignal
% index_compute_raw = knnsearch(target,source,'k',1);
% 
% index_gt = 1:size(target,1);
% 
% comp_error = index_compute_raw - index_gt';
% 
% index_count = find(comp_error==0);
% 
% rate = length(index_count);
% 
% fprintf('original accurate rate is %d percent\n',100*rate/size(index_compute_raw,1));



end

