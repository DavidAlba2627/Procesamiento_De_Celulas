close all;
clear variables;
clc;

% Leyendo la imagen y separando por canal de color
cd ..; cd cell_images;
Im1 = im2double(imread('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(b).tif'));
cd ..; cd scripts;
figure();imshow(Im1); title('Imagen original','Fontsize',15);
Im1_C = Im1(:,:,3); 
figure();imshow(Im1_C); title('Imagen de las celulas','Fontsize',15);
area_celula = 400;

aux = Im1_C > 0.1;
Im1_C = aux.*Im1_C;
figure();imshow(Im1_C); title('Umbral global','Fontsize',15);

% Realizando el umbral local
g = localthresh(Im1_C, ones(9), 0, 0.9);
figure(); imshow(g, []);
% Rellenando los agujeros
Celulas_Rellenas = imfill(g,'holes');
figure(); imshow(Celulas_Rellenas, []); title('Imagen despues de rellenar agujeros');
% Realizando la apertura
B = strel('disk', 5);
Celulas_Apertura = imopen(Celulas_Rellenas , B);
figure(); imshow(Celulas_Apertura, []); title('Imagen despues de realizar la apertura');
% Calculando propiedades
[L, n, stats, Tabla] = Calcular_Propiedades(Celulas_Apertura);
Area = stats.Area;
% Intensidad media
Inten_media_nucleo = cell2mat(struct2cell(regionprops(L, Im1_C, 'MeanIntensity')));
Condiciones_Intensidad = [Inten_media_nucleo'] < 0.06;
cc = bwconncomp(Celulas_Apertura); 
idx = find(Condiciones_Intensidad); 
BW = ismember(labelmatrix(cc),idx);  
Celulas_Inten_Correcta = Celulas_Apertura - BW;
figure(); imshow(Celulas_Inten_Correcta) ; title('Celulas con intensidades correctas','Fontsize',15);
% Calculando propiedades
[L1, n1, stats1, Tabla1] = Calcular_Propiedades(Celulas_Inten_Correcta);

% Imagen de las celulas enumeradas
figure();imshow(Celulas_Inten_Correcta); title('Celulas enumeradas depuradas en inten','Fontsize',15);
hold on;
for k=1:n
    [r,c] = find(L1==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end

% Depurando por area pequeña
Condicion_Area = [stats1.Area] < 150;
Celulas_Area_Correcta = Eliminar_Nucleos(Celulas_Inten_Correcta, Condicion_Area);
[L2, n2, stats2, Tabla2] = Calcular_Propiedades(Celulas_Area_Correcta);
figure(); imshow(Celulas_Area_Correcta) ; title('Celulas con area correcta','Fontsize',15);
figure();imshow(Celulas_Area_Correcta); title('Celulas enumeradas depuaradas en area','Fontsize',15);
hold on;
for k=1:n2
    [r,c] = find(L2==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end

% Dividiendo las celulas
Celulas_Divididas = Celulas_Area_Correcta;
% Repitiendo varias veces la separacion de celulas
% for i = 1:3
%     % Obteniendo las propiedades de los nucleos recuperados
%     [L3, n3, stats3, Tabla3] = Calcular_Propiedades(Celulas_Divididas);
%     % Separando nucleos
%     Condiciones_Separacion2 = [stats3.MajorAxisLength./stats3.MinorAxisLength] > 1.5 & [stats3.Area] > area_celula*1.5;
%     Celulas_Divididas = Separar_Nucleos(Celulas_Divididas, Condiciones_Separacion2);
%     figure();imshow(Celulas_Divididas); title('Segmentacion de las celulas','Fontsize',15);
% end

[L3, n3, stats3, Tabla3] = Calcular_Propiedades(Celulas_Divididas);
% Separando nucleos
Condiciones_Separacion2 = [stats3.MajorAxisLength./stats3.MinorAxisLength] > 1.4 & [stats3.Area] > area_celula*1.5 & [stats3.Area] < area_celula*4.0 & [stats3.Circularity] > 0.4;
Celulas_Divididas = Separar_Nucleos(Celulas_Divididas, Condiciones_Separacion2);
figure();imshow(Celulas_Divididas);title('Division de celulas','Fontsize',15);



figure();imshow(Celulas_Divididas); title('Segmentacion de celulas','Fontsize',15);
% Separando las celulas en duda
[L4, n4, stats4, Tabla4] = Calcular_Propiedades(Celulas_Divididas);
Condiciones_Marcado_Rojo = [stats4.Area] > area_celula*3;
Celulas_Correctas_Seg = Eliminar_Nucleos(Celulas_Divididas, Condiciones_Marcado_Rojo);
Celulas_Duda = Celulas_Divididas - Celulas_Correctas_Seg;
figure();imshow(Celulas_Correctas_Seg);title('Celulas correctas segmentadas','Fontsize',15);
figure();imshow(Celulas_Duda);title('Celulas en duda','Fontsize',15);


% Imagen de las celulas enumeradas
figure();imshow(Celulas_Divididas); title('Celulas enumeradas despues de la division','Fontsize',15);
hold on;
for k=1:n4
    [r,c] = find(L4==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end

% Imagen en azul
[r,c] = size(Celulas_Divididas);
etiquetas = zeros(r,c,3);
etiquetas(:,:,3) =  Celulas_Divididas;
figure();imshow(etiquetas);


% Segmentacion de canales de color
% Aplicando el métod de Otsu
[ImR_seg, ImR_gris] = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(r).tif',1,1);
[ImG_seg, ImG_gris] = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(g).tif',2,1);
[ImIR_seg, ImIR_gris] = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(IR).tif',0,1);

[Im_Class_R, Tabla_Propiedades_R] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, ImR_seg,  Im1_C, ImR_gris);
title('Clasificacion sobre el canal ROJO','Fontsize',15);
[Im_Class_G, Tabla_Propiedades_G] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, ImG_seg,  Im1_C, ImG_gris);
title('Clasificacion sobre el canal VERDE','Fontsize',15);
[Im_Class_IR, Tabla_Propiedades_IR] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, ImIR_seg,  Im1_C, ImIR_gris);
title('Clasificacion sobre el canal INFRARROJO','Fontsize',15);
% Inten_media_nucleo = cell2mat(struct2cell(regionprops(L, Im1_C, 'MeanIntensity')));
% Condiciones_Intensidad = [Inten_media_nucleo'] < 0.06;
% cc = bwconncomp(Celulas_Apertura); 
% idx = find(Condiciones_Intensidad); 
% BW = ismember(labelmatrix(cc),idx);  
% Celulas_Correctas = Celulas_Apertura - BW;
% figure(); imshow(Celulas_Correctas) ; title('Celulas corrctas','Fontsize',15);

ID_G = Tabla_Propiedades_G(:, 'Estado');
ID_G = ID_G{:,:};
ID_IR = Tabla_Propiedades_IR(:, 'Estado');
ID_IR = ID_IR{:,:};

% Etiquetas de los canales verde e infrarrojo
G_Pos = ID_G == 1;
G_Neg = ID_G == 0;
IR_Pos = ID_IR == 1;
IR_Neg = ID_IR == 0;
% Ambos positivos
GP_IRP = G_Pos & IR_Pos;
% Verde positivo e infrarrojo negativo
GP_IRN = G_Pos & IR_Neg;
% Verde negativo e infrarrojo positivo
GN_IRP = G_Neg & IR_Pos;
% Ambos negativos
GN_IRN = G_Neg & IR_Neg;

Im_GP_IRP = zeros(r,c);
Im_GP_IRN = zeros(r,c);
Im_GN_IRP = zeros(r,c);
Im_GN_IRN = zeros(r,c);

for i=1:n4
    Celula = L4==i;
    Im_GP_IRP = Im_GP_IRP + Celula*GP_IRP(i); 
    Im_GP_IRN = Im_GP_IRN + Celula*GP_IRN(i); 
    Im_GN_IRP = Im_GN_IRP + Celula*GN_IRP(i); 
    Im_GN_IRN = Im_GN_IRN + Celula*GN_IRN(i);
end


figure();imshow(Im_GP_IRP); title('Verde positivo e infrarrojo positivo','Fontsize',15);
figure();imshow(Im_GP_IRN); title('Verde positivo e infrarrojo negativo','Fontsize',15);
figure();imshow(Im_GN_IRP); title('Verde negativo e infrarrojo positivo','Fontsize',15);
figure();imshow(Im_GN_IRN); title('Verde negativo e infrarrojo negativo','Fontsize',15);
% figure();imshow(Im_GP_IRP + Im_GP_IRN + Im_GN_IRP + Im_GN_IRN); title('Verde negativo e infrarrojo negativo','Fontsize',15);
