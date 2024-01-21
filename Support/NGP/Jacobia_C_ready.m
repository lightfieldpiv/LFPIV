function Matrix = Jacobia_C_ready ( x,control_point,source_obj,...
                                     weight,weightTrans,...
                                     coeff_rigid ,coeff_smooth,alpha_point )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% m_size = size(control_point,1)*6 + size(control_point,1)^2*3 + size(source_obj,1)*3 + size(source_obj,1);
m_size = size(control_point,1)*6 + size(control_point,1)^2*3 + size(source_obj,1)*3;

% rigid part, smooth part, fit part respectively

variable_num = 12;

% Matrix = zeros(m_size,variable_num*size(control_point,1));
Matrix = sparse(m_size,variable_num*size(control_point,1));

x_shape = reshape(x,12,size(control_point,1));

a1 = x_shape(1:3,:);
a2 = x_shape(4:6,:);
a3 = x_shape(7:9,:);
% b  = x_shape(10:12,:);

%---------- seperate to x1,x2,x3,,,,,,
x1 = a1(1,:); x4 = a2(1,:); x7 = a3(1,:); %x10 = b(1,:);
x2 = a1(2,:); x5 = a2(2,:); x8 = a3(2,:); %x11 = b(2,:);
x3 = a1(3,:); x6 = a2(3,:); x9 = a3(3,:); %x12 = b(3,:);

%% ----- rigid part  1----------
   
%----- rigid part  1----------
seq = 1;
row_index = 1; 
for i = size(control_point,1)*(seq-1)+1 : size(control_point,1)*seq
    % first 6 point has jabobia
    index1 = (i-1)*variable_num +1;
    %index2 = (i-1)*variable_num +2;
    %index3 = (i-1)*variable_num +3;
    %index4 = (i-1)*variable_num +4;
    %index5 = (i-1)*variable_num +5;
    index6 = (i-1)*variable_num +6;
    %-------
    Matrix(i, index1:index6) = [x4(i),x5(i),x6(i),x1(i),x2(i),x3(i)];
    %Matrix(i, index2) = x5(i);
    %Matrix(i, index3) = x6(i);
    
    %Matrix(i, index4) = x1(i);
    %Matrix(i, index5) = x2(i);
    %Matrix(i, index6) = x3(i);
    
    row_index = row_index + 1; 
end
%----------------------
seq = seq + 1;
%------ part 2

% for i = size(control_point,1)*(seq-1)+1 : size(control_point,1)*seq
for i = 1 : size(control_point,1)
    
    %row_index = size(control_point,1)*(seq-1)+i;
    % first 6 point has jabobia
    index1 = (i-1)*variable_num +1;
    %index2 = (i-1)*variable_num +2;
    index3 = (i-1)*variable_num +3;
    index7 = (i-1)*variable_num +7;
    %index8 = (i-1)*variable_num +8;
    index9 = (i-1)*variable_num +9;
    %-------
    Matrix(row_index, index1:index3) = [x7(i),x8(i),x9(i)];
    %Matrix(row_index, index2) = x8(i);
   % Matrix(row_index, index3) = x9(i);
    
    Matrix(row_index, index7:index9) = [x1(i),x2(i),x3(i)];
%     Matrix(row_index, index8) = x2(i);
%     Matrix(row_index, index9) = x3(i);
    
    row_index = row_index + 1; 
end
seq = seq + 1;
%------ part 3

for i = 1 : size(control_point,1)
    
    %row_index = size(control_point,1)*(seq-1)+i;
    % first 6 point has jabobia
    index4 = (i-1)*variable_num +4;
    %index5 = (i-1)*variable_num +5;
    %index6 = (i-1)*variable_num +6;
    %index7 = (i-1)*variable_num +7;
    %index8 = (i-1)*variable_num +8;
    index9 = (i-1)*variable_num +9;
    %-------
    Matrix(row_index, index4:index9) = [x7(i),x8(i),x9(i),x4(i),x5(i),x6(i)];
%     Matrix(row_index, index5) = x8(i);
%     Matrix(row_index, index6) = x9(i);
%     
%     Matrix(row_index, index7) = x4(i);
%     Matrix(row_index, index8) = x5(i);
%     Matrix(row_index, index9) = x6(i);
    
    row_index = row_index + 1; 
end
seq = seq + 1;
% ----- rigid part  2 ----------
for i = 1 : size(control_point,1)
    
    %row_index = size(control_point,1)*(seq-1)+i;
    % first 6 point has jabobia
    index1 = (i-1)*variable_num +1;
    %index2 = (i-1)*variable_num +2;
    index3 = (i-1)*variable_num +3;
    
    %-------
    Matrix(row_index, index1:index3) = [-2 * x1(i),-2 * x2(i),-2 * x3(i)];
   % Matrix(row_index, index2) = -2 * x2(i);
    %Matrix(row_index, index3) = -2 * x3(i);
    
     row_index = row_index + 1; 
end
seq = seq + 1;

%---- rigid part 2.2
for i = 1 : size(control_point,1)
    
   % row_index = size(control_point,1)*(seq-1)+i;
    % first 6 point has jabobia
    index4 = (i-1)*variable_num + 4;
    index5 = (i-1)*variable_num + 5;
    index6 = (i-1)*variable_num + 6;
    
    %-------
    Matrix(row_index, index4) = -2 * x4(i);
    Matrix(row_index, index5) = -2 * x5(i);
    Matrix(row_index, index6) = -2 * x6(i);
    
     row_index = row_index + 1;
end
seq = seq + 1;

%---- rigid part 2.3
for i = 1 : size(control_point,1)
    
    %row_index = size(control_point,1)*(seq-1)+i;
    % first 6 point has jabobia
    index7 = (i-1)*variable_num + 7;
    index8 = (i-1)*variable_num + 8;
    index9 = (i-1)*variable_num + 9;
    
    %-------
    Matrix(row_index, index7) = -2 * x7(i);
    Matrix(row_index, index8) = -2 * x8(i);
    Matrix(row_index, index9) = -2 * x9(i);
    
    row_index = row_index + 1;
end

row_index = row_index-1;
Matrix(1:row_index,:) =  sqrt(coeff_rigid) * Matrix(1:row_index,:);
seq = seq + 1;
%% E smooth term
control_trans = control_point';
control_point_num = size(control_point,1);

row_smooth = size(control_point,1)*(seq-1)+1;
smooth_start = size(control_point,1)*(seq-1)+1;
smooth_end  = size(control_point,1)*(seq-1)+control_point_num^2*3;

coeff =  sqrt(coeff_smooth);
for i = 1 : control_point_num
        index1i  = (i-1)*variable_num + 1;
        index2i  = (i-1)*variable_num + 2;
        index3i  = (i-1)*variable_num + 3;
        index4i  = (i-1)*variable_num + 4;
        index5i  = (i-1)*variable_num + 5;
        index6i  = (i-1)*variable_num + 6;
        index7i  = (i-1)*variable_num + 7;
        index8i  = (i-1)*variable_num + 8;
        index9i  = (i-1)*variable_num + 9;
        index10i = (i-1)*variable_num + 10;
        index11i = (i-1)*variable_num + 11;
        index12i = (i-1)*variable_num + 12;
        %--------
    for j = 1 : control_point_num
       
        x_i  = control_trans(:,i);
        x_j  = control_trans(:,j);
        D    = x_j - x_i;
        Dx   = D(1,:);
        Dy   = D(2,:);
        Dz   = D(3,:);
        eff  = sqrt(weight(i,j));
        %-------
        index10j = (j-1)*variable_num + 10;
        index11j = (j-1)*variable_num + 11;
        index12j = (j-1)*variable_num + 12;
        %--------
         for k = 1 : 3
        %---- ever iter plus one
            if(k == 1) 
                 Matrix(row_smooth, index1i)  = coeff * Dx * eff;
                 Matrix(row_smooth, index4i)  = coeff * Dy * eff;
                 Matrix(row_smooth, index7i)  = coeff * Dz * eff;
                 Matrix(row_smooth, index10i) = coeff * 1 * eff;
                 Matrix(row_smooth, index10j) = coeff *-1 * eff;
            elseif(k==2)
                 Matrix(row_smooth, index2i)  =  coeff *Dx * eff;
                 Matrix(row_smooth, index5i)  =  coeff *Dy * eff;
                 Matrix(row_smooth, index8i)  =  coeff *Dz * eff;
                 Matrix(row_smooth, index11i) =  coeff *1 * eff;
                 Matrix(row_smooth, index11j) = coeff *-1 * eff;
            elseif(k==3)
                 Matrix(row_smooth, index3i)  =  coeff *Dx * eff;
                 Matrix(row_smooth, index6i)  =  coeff *Dy * eff;
                 Matrix(row_smooth, index9i)  =  coeff *Dz * eff;
                 Matrix(row_smooth, index12i) =  coeff *1 * eff;
                 Matrix(row_smooth, index12j) = coeff *-1 * eff;
            end
%          Matrix(row_smooth, index2i) = Matrix(row_smooth, index1i);
%          Matrix(row_smooth, index3i) = Matrix(row_smooth, index1i);
%          
%          Matrix(row_smooth, index4i) = Dy * eff;
%          Matrix(row_smooth, index5i) = Matrix(row_smooth, index4i);
%          Matrix(row_smooth, index6i) = Matrix(row_smooth, index4i);
%          
%          Matrix(row_smooth, index7i) =  Dz * eff;
%          Matrix(row_smooth, index8i) =  Matrix(row_smooth, index7i);
%          Matrix(row_smooth, index9i) =  Matrix(row_smooth, index7i);
%          
%          Matrix(row_smooth, index10i) = 1 * eff;
%          Matrix(row_smooth, index11i) = Matrix(row_smooth, index10i);
%          Matrix(row_smooth, index12i) = Matrix(row_smooth, index10i);
%          %-----
%          Matrix(row_smooth, index10j) = -1 * eff;
%          Matrix(row_smooth, index11j) = Matrix(row_smooth, index10j);
%          Matrix(row_smooth, index12j) = Matrix(row_smooth, index10j);
         
         row_smooth = row_smooth + 1;
        end
    end
end

% Matrix(smooth_start:smooth_end,:) =  sqrt(coeff_smooth) * Matrix(smooth_start:smooth_end,:);

%% E fit_1 term
% alpha_point = 0.1;
% alpha_plane = 1.0;
for i = 1 : size(source_obj,1)
    
%         eff1 = 0; eff4 = 0; eff7 = 0; eff10 = 0;
%         eff2 = 0; eff5 = 0; eff8 = 0; eff11 = 0;
%         eff3 = 0; eff6 = 0; eff9 = 0; eff12 = 0;
        %--------------
         alpha = sqrt(alpha_point);
    for k = 1 : 3 
        for j = 1 : control_point_num
                D = source_obj(i,:)' - control_point(j,:)';
                Dx   = D(1,:); Dy = D(2,:); Dz   = D(3,:);
                %-----------------------------
                index1j  = (j-1)*variable_num + 1;
                index2j  = (j-1)*variable_num + 2;
                index3j  = (j-1)*variable_num + 3;
                index4j  = (j-1)*variable_num + 4;
                index5j  = (j-1)*variable_num + 5;
                index6j  = (j-1)*variable_num + 6;
                index7j  = (j-1)*variable_num + 7;
                index8j  = (j-1)*variable_num + 8;
                index9j  = (j-1)*variable_num + 9;
                index10j = (j-1)*variable_num + 10;
                index11j = (j-1)*variable_num + 11;
                index12j = (j-1)*variable_num + 12;
                %-----
                if(k==1)
                    Matrix(row_smooth, index1j)  = alpha*weightTrans(i,j) * Dx;
                    Matrix(row_smooth, index4j)  = alpha*weightTrans(i,j) * Dy;
                    Matrix(row_smooth, index7j)  = alpha*weightTrans(i,j) * Dz;
                    Matrix(row_smooth, index10j) = alpha*weightTrans(i,j) * 1;
                elseif(k==2)
                    Matrix(row_smooth, index2j)  = alpha*weightTrans(i,j) * Dx;
                    Matrix(row_smooth, index5j)  = alpha*weightTrans(i,j) * Dy;
                    Matrix(row_smooth, index8j)  = alpha*weightTrans(i,j) * Dz;
                    Matrix(row_smooth, index11j) = alpha*weightTrans(i,j) * 1;
                elseif(k==3)
                    Matrix(row_smooth, index3j)  = alpha*weightTrans(i,j) * Dx;
                    Matrix(row_smooth, index6j)  = alpha*weightTrans(i,j) * Dy;
                    Matrix(row_smooth, index9j)  = alpha*weightTrans(i,j) * Dz;
                    Matrix(row_smooth, index12j) = alpha*weightTrans(i,j) * 1;
                end
%                  Matrix(row_smooth, index1j)  = alpha*weightTrans(i,j) * Dx;
%                  Matrix(row_smooth, index2j)  = alpha*weightTrans(i,j) * Dx;
%                  Matrix(row_smooth, index3j)  = alpha*weightTrans(i,j) * Dx;
% 
%                  Matrix(row_smooth, index4j)  = alpha*weightTrans(i,j) * Dy;
%                  Matrix(row_smooth, index5j)  = alpha*weightTrans(i,j) * Dy;
%                  Matrix(row_smooth, index6j)  = alpha*weightTrans(i,j) * Dy;
% 
%                  Matrix(row_smooth, index7j)  = alpha*weightTrans(i,j) * Dz;
%                  Matrix(row_smooth, index8j)  = alpha*weightTrans(i,j) * Dz;
%                  Matrix(row_smooth, index9j)  = alpha*weightTrans(i,j) * Dz;
% 
%                  Matrix(row_smooth, index10j) = alpha*weightTrans(i,j) * 1;
%                  Matrix(row_smooth, index11j) = alpha*weightTrans(i,j) * 1;
%                  Matrix(row_smooth, index12j) = alpha*weightTrans(i,j) * 1;

        end
             row_smooth = row_smooth + 1;  
    end
end


end

