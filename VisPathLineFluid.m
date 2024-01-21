function VisPathLineFluid(uv, param, folder)
H_size = param.H_size+1;
W_size = param.W_size+1;
D_size = param.D_size+1;

[x, y, z] = meshgrid(W_size:-1:1,H_size:-1:1,D_size:-1:1);
% [x, y, z] = meshgrid(1:1:W_size,1:1:H_size,1:1:D_size);

coefficient = 8;

for i = 1 : size(uv,1)
%     fid = fopen('Bvec.raw','w'); fwrite(fid,uv{i},'float'); fclose(fid);
    uv{i}(:,:,:,1) =  coefficient*uv{i}(:,:,:,1);
    uv{i}(:,:,:,2) =  coefficient*uv{i}(:,:,:,2);
    uv{i}(:,:,:,3) =  coefficient*uv{i}(:,:,:,3);

    u = uv{i}(:,:,:,1);
    v = uv{i}(:,:,:,2);
    w = uv{i}(:,:,:,3);
    
    mag_field = Magnitude(uv{i});
    
    c_name = sprintf('flow_%d.vtk',i);
    
    vtkwrite([folder,c_name], 'structured_grid', x, y, z, ... 
    'vectors', 'vector_field', u, v, w,'scalars','mag_field',mag_field);
    
%     B = cat(4,u,v,w); 
%     B = permute(B,[4 1 2 3]); 
%     
%     fid = fopen('Bvec.raw','w'); fwrite(fid,B,'float'); fclose(fid);
end