function mag_field = Magnitude(uv)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

u = uv(:,:,:,1);
v = uv(:,:,:,2);
w = uv(:,:,:,3);

mag_field = zeros(size(u));

for i = 1 : size(u,1)
    for j = 1 : size(u,2)
     for k = 1 : size(u,3)
            mag_field(i,j,k) =  sqrt(u(i,j,k)^2 +  v(i,j,k)^2  + w(i,j,k)^2 );
     end
    end
end

