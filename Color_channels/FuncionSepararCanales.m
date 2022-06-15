function [] = FuncionSepararCanales(File_name)
% Leemos imagen usando el script de la compania Zeiss
% File_name = 'Taiep2_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)';
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
end

