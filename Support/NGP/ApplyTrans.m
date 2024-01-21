function result = ApplyTrans(source_obj,a1,a2,a3,b,control_point,weightTrans)
%UNTITLED7 Summary of this function goes here

source_length  = size(source_obj,1);
control_length = size(control_point,1);

%%
source_trans = source_obj';
control_trans = control_point';
result = zeros(3, source_length);
for i = 1 : control_length
    
    A_i    = [a1(:,i),a2(:,i),a3(:,i)];
    x_i    = repmat(control_trans(:,i),1,source_length);
    b_i    = repmat(b(:,i),1,source_length);
    value  = A_i * (source_trans - x_i) + x_i + b_i;
    value  = repmat(weightTrans(:,i)',3,1).*value; % multiple by weight
    result = result + value;
   
end

result = result';

end

