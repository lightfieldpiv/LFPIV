function s_resize = scalePoints(s_pos,min_x,max_x,...
                                      min_y,max_y,...
                                      min_z,max_z,...
                                      H_size,...
                                      W_size,...
                                      D_size)
s_resize = zeros(size(s_pos));
for i = 1 : size(s_pos,1)
    s_resize(i,1) = (s_pos(i,1) - min_x)/(max_x - min_x) * W_size + 1;
    s_resize(i,2) = (s_pos(i,2) - min_y)/(max_y - min_y) * H_size + 1;
    s_resize(i,3) = (s_pos(i,3) - min_z)/(max_z - min_z) * D_size + 1;

end


end