function [Im_Seg, Im1_gris] = Seg_canales(nombre,n, verbose)
% Segmenta las celulas con el método de otsu
% n=1 todas las imagenes, n=2 solo la imagen segmentada n=0 no muestra ninnguna imagen 
cd ..; cd cell_images;
Im1 = im2double(imread(nombre));
cd ..; cd scripts;

% if canal == 0
%     Im1_gris = Im1;
% else
%     Im1_gris = Im1(:,:,canal);
% end

Im1_gris = Im1;
[h,b] = imhist(Im1_gris, 256);
[T,SM] = otsuthresh(h);
Im_Seg = im2bw(Im1_gris,T);
if verbose

    if n == 1
    figure();imshow(Im1); title(' Canal Infrarrojo','Fontsize',15);
    figure();imshow(Im_Seg);title(' Segmentacion del canal Infrarrojo','Fontsize',15);
    elseif n == 2
    figure();imshow(Im_Seg);title(' Segmentacion del canal Infrarrojo','Fontsize',15);
    end
        
end
end

