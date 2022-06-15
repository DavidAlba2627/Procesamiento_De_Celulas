% Probando funcion Zeiss

clc; 
%close all; 
clear variables;

% Leemos imagen usando el script de la compania Zeiss
File_name = 'Taiep2_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)';
out = ReadImage6D(strcat('Images/', File_name, '.czi'), false, 0);
metadata = out{2};
image6d = out{1};

%Nombre de marcadores
Dyes = metadata.Dyes;
%Longitudes de emision
WLEm = metadata.WLEm;
%Encontrando indice de DAPI
Ind_Dapi = find(strcmp(Dyes, 'DAPI'));

%Encontrando escalas
Escala = [metadata.ScaleX, metadata.ScaleY];

%Encontrando los demas indices
Indices_Num = zeros(4,1);
Indices_Nombre = string(zeros(4,1));
Indices_Nombre_Sin_DAPI = string(zeros(4,1));
figure();
for i=1:size(Dyes,2)
    if Dyes{i} == "DAPI"
        borrar = i;
    end
    Indices_Nombre(i) = Dyes{i};
    Indices_Nombre_Sin_DAPI(i) = Dyes{i};
    Indices_Num(i) = i; 
    %Separamos la imagen en canales de color
    img_ch = image6d(1,1,1,i,:,:);
    img_ch = squeeze(img_ch);
    % Reescalamos los datos de las intensidades
    img_ch = img_ch - min(img_ch(:)); img_ch = img_ch / max(img_ch(:));
    subplot(1,size(Dyes,2),i); imshow(img_ch, []); title(Indices_Nombre(i));
    cd ..; cd cell_images;
    imwrite(img_ch, strcat(File_name, Dyes{i}, '.tif'));
    cd ..; cd Color_channels;
        
    
end

Indices_Nombre_Sin_DAPI(borrar,:) = [];



cd ..; cd cell_images;
save(strcat(File_name,'Canales.mat'),'Indices_Nombre', 'Indices_Nombre_Sin_DAPI', 'Escala');
%Separamos la imagen en canales de color
% img_ch1 = image6d(1,1,1,2,:,:);
% img_ch2 = image6d(1,1,1,1,:,:);
% img_ch3 = image6d(1,1,1,4,:,:);
% img_ch4 = image6d(1,1,1,3,:,:);

% img_ch1 = image6d(1,1,1,1,:,:);
% img_ch2 = image6d(1,1,1,2,:,:);
% img_ch3 = image6d(1,1,1,3,:,:);
% img_ch4 = image6d(1,1,1,4,:,:);
% img_ch1 = squeeze(img_ch1);
% img_ch2 = squeeze(img_ch2);
% img_ch3 = squeeze(img_ch3);
% img_ch4 = squeeze(img_ch4);

% figure()
% subplot(1,4,1); imagesc(img_ch1); axis equal; axis tight;
% subplot(1,4,2); imagesc(img_ch2); axis equal; axis tight;
% subplot(1,4,3); imagesc(img_ch3); axis equal; axis tight;
% subplot(1,4,4); imagesc(img_ch4); axis equal; axis tight;

% % Reescalamos los datos de las intensidades
% img_ch1 = img_ch1 - min(img_ch1(:)); img_ch1 = img_ch1 / max(img_ch1(:));
% img_ch2 = img_ch2 - min(img_ch2(:)); img_ch2 = img_ch2 / max(img_ch2(:));
% img_ch3 = img_ch3 - min(img_ch3(:)); img_ch3 = img_ch3 / max(img_ch3(:));
% img_ch4 = img_ch4 - min(img_ch4(:)); img_ch4 = img_ch4 / max(img_ch4(:));
% 
% 
% 
% figure()
% %Subplot de los canales de color
% %Primero en escala de grises
% subplot(2,4,1); imshow(img_ch1, []); title('canal 1 - verde');
% subplot(2,4,2); imshow(img_ch2, []); title('canal 2 - IR');
% subplot(2,4,3); imshow(img_ch3, []); title('canal 3 - Azul');
% subplot(2,4,4); imshow(img_ch4, []); title('canal 4 - Rojo');

% % A las escalas de grises agregamos colores tal como indicaron los expertos
% [m, n, r] = size(img_ch1);
% % Canal 1 asignamos verde
% img_g=zeros(m,n,3);
% img_g(:,:,2)=img_ch1;
% subplot(2,4,5),imshow(img_g, []); title('canal 1 - verde');
% % Canal 2 lo dejamos en gris como convencion para el IR 
% imh_IR = img_ch2;
% subplot(2,4,6),imshow(imh_IR, []); title('canal 2 - IR');
% % Canal 3 asignamos azul
% img_b=zeros(m,n,3);
% img_b(:,:,3)=img_ch3;
% subplot(2,4,7),imshow(img_b, []); title('canal 3 - Azul');
% %Canal 4 asignamos rojo
% img_r=zeros(m,n,3);
% img_r(:,:,1)=img_ch4;
% subplot(2,4,8),imshow(img_r, []); title('canal 4 - Rojo');


% % Nos cambioamos de directorio para guardar las imagenes
% cd ..; cd cell_images;
% % Guardamos las imagenes de los canales de color
% imwrite(img_ch1, strcat(File_name, '(g) .tif'))
% imwrite(img_ch2,strcat(File_name, '(gray).tif'))
% imwrite(img_ch3, strcat(File_name, '(b).tif'))
% imwrite(img_ch4, strcat(File_name, '(r).tif'))