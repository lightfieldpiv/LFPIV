function cnt = detectParticleLoc(img,thres)
% img_name = 'IMG_0223__Decoded_Thumb.png';
a =  img;
% thres = 30;
% a = rgb2gray(a);
% imshow(a,[])
% a = a*255;
% a = 255-a;
a = a(:,:,2);
a = a*255;

% colormap('gray');figure;imshow(a,[]);
b = bpass(a,1,10);
% colormap('gray');figure; imshow(b,[]);
% max(max(b))
pk = pkfnd(b,thres,11);

cnt = cntrd(b,pk,15);
hist(mod(cnt(:,1),1),20);

radi = 3*ones(size(cnt,1),1);
imshow(b,[]);
h = viscircles(cnt(:,1:2),radi)

end
% h = viscircles(centers,radii);
% ob_img = imfilter(img_current,fspecial('Gaussian',8, 1.2));
% whos a
% imshow(a);
% a = 255-a;
% colormap('gray'), imagesc(a);