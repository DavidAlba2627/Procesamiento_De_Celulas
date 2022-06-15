close all;
clear variables;
clc;

% Segundo commit xd
% Leyendo la imagen y separando por canal de color
cd ..; cd cell_images;
Im1 = im2double(imread('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(b).tif'));
cd ..; cd scripts;

figure();imshow(Im1); title(' Imagen original canal Azul','Fontsize',15);
Im1_C = Im1(:,:,3);


figure();imshow(Im1_C); title(' Imagen de las celulas canal Azul','Fontsize',15);

% Aplicando el métod de Otsu
ImR_seg = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(r).tif',1,1);
ImG_seg = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(g).tif',2,1);
ImIR_seg = Seg_canales('SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)(IR).tif',0,1);
    
% Segmentando los nucleos
Celulas = Im1_C > 15/255;
B1 = strel('disk', 5);  
Im1_C_Seg = imopen(Celulas, B1); 
Im1_C_Seg = imfill(Im1_C_Seg,'holes');
figure();imshow(Im1_C_Seg);title('1. Segmentacion de las celulas','Fontsize',15);




% Obteniendo las propiedades de los nucleos
[L1, n1, stats1, Tabla1] = Calcular_Propiedades(Im1_C_Seg);
figure();imshow(Im1_C); title('Celulas enumeradas antes de las condiciones','Fontsize',15);
hold on;
for k=1:n1
    [r,c] = find(L1==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end
% Eliminando nucleos que no nos importan (Grandes y no circulares)
area_condicion = 450;
Condiciones = [stats1.Area] < area_condicion/3 | [stats1.Area] > area_condicion*3 | [stats1.Circularity] < 0.3;
Celulas_Depuradas = Eliminar_Nucleos(Im1_C_Seg, Condiciones);
figure();imshow(Celulas_Depuradas); title('Celulas buenas','Fontsize',15);

% Obteniendo las propiedades de los nucleos depurados
[L2, n2, stats2, Tabla2] = Calcular_Propiedades(Celulas_Depuradas);
% Separando nucleos
Condiciones_Separacion = [stats2.MajorAxisLength./stats2.MinorAxisLength] > 1.5 & [stats2.Area] > area_condicion*1.5;
Celulas_Separadas = Separar_Nucleos(Celulas_Depuradas, Condiciones_Separacion);
figure();imshow(Celulas_Separadas); title('Segmentacion de celulas','Fontsize',15);

% Celulas ovaladas
[L3, n3, stats3, Tabla3] = Calcular_Propiedades(Celulas_Separadas);
Condiciones_Ovaladas = [stats3.Circularity] < 0.5;
Celulas_No_Ovaladas = Eliminar_Nucleos(Celulas_Separadas, Condiciones_Ovaladas);
Celulas_Ovaladas = Celulas_Separadas - Celulas_No_Ovaladas;
% Crecimiento de regiones
[r, c] = size(Celulas_Ovaladas);
Matriz = ones(r,c);
Umbral = 1.5*Im1_C;
[g, NR, SI, TI] = regiongrow(Im1_C, Celulas_Ovaladas, Umbral);
aux = g>=1;
B2 = strel('disk', 10);
Celulas_Crecidas = imclose(aux, B2);
%figure(); imshow(Celulas_Crecidas); title('Imagen de crecimiento de regiones','Fontsize',15);
% Eliminando celulas que crecieron mucho
[L4, n4, stats4, Tabla4] = Calcular_Propiedades(Celulas_Crecidas);
Condiciones_Crecidas = [stats4.Area] < area_condicion/3 | [stats4.Area] > area_condicion*2 | [stats4.Circularity] < 0.3;
Celulas_No_Crecidas = Eliminar_Nucleos(Celulas_Crecidas, Condiciones_Crecidas);
Segmentacion_Celulas = Celulas_No_Ovaladas + Celulas_No_Crecidas; 
figure(); imshow(Segmentacion_Celulas); title('Segmemtacion de celulas','Fontsize',15);

% ------ Proteina (antes farmaco) ---------
%Im1_F_Seg = Im1_F > 0.14;
% figure();imshow(Im1_F_Seg, []); title('Segmentacion de la proteina','Fontsize',15);

% Tabla de propiedades de las celulas
[L5, n5, stats5, Tabla5] = Calcular_Propiedades(Segmentacion_Celulas);
Num5 = [1:n5]';
area5 = stats5.Area;
eje_Mayor5 = stats5.MajorAxisLength;
eje_Menor5 = stats5.MinorAxisLength;
orientacion5 = stats5.Orientation;
Inten_media_nucleo = cell2mat(struct2cell(regionprops(L5, Im1_C, 'MeanIntensity')));
%Inten_media_proteina = cell2mat(struct2cell(regionprops(L5, Im1_F, 'MeanIntensity')));

% % Calculando positivos/negativos
% figure();
% h_Im = histogram(L5,n5+1);
% h_Im_val = h_Im.Values;
% h_Im_val = h_Im_val(2:end);
% p_o_n = zeros(length(n5));
% [r,c] = size(L4);
% Im_rgb = zeros(r,c,3);
% porcentaje_a = zeros(n5,1);
% a_protenia = zeros(n5,1);
% for i = 1:n5
%      Celula = L5 == i;
%      Area_C = h_Im_val(i);
%      Superposicion = Celula.*Im1_F_Seg;
%      Area_P = sum(Superposicion(:));
%      
%      if Area_P < .3*Area_C
%        p_o_n(i) = 0;
%        Im_rgb(:,:,3) = Im_rgb(:,:,3) + Celula;
%      else
%        p_o_n(i) = 1;
%        Im_rgb(:,:,2) = Im_rgb(:,:,2) + Celula;
%      end
%      porcentaje_a(i,1) = (Area_P*100)/Area_C;
%      a_protenia(i,1) = Area_P;
% end
% 
% figure();imshow(Im_rgb);title('Positivos y Negativos','Fontsize',15);
% Tabla_Propiedades = table(Num5, area5, a_protenia, porcentaje_a, eje_Mayor5, eje_Menor5, orientacion5, Inten_media_nucleo', Inten_media_proteina', p_o_n');
% Tabla_Propiedades.Properties.VariableNames = {'Elemento','Area','Area proteina dentro de la celuala','Porcentaje', 'Eje mayor', 'Eje menor', 'Orientacion', 'Inten_media nucleo', 'Inten_media_proteina','Estado'};


% Imagen de las celulas enumeradas
figure();imshow(Segmentacion_Celulas); title('Celulas enumeradas','Fontsize',15);
hold on;
for k=1:n5
    [r,c] = find(L5==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end
figure();imshow(Im1_C); title('Celulas enumeradas','Fontsize',15);
hold on;
for k=1:n5
    [r,c] = find(L5==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'r');
end


% Umbralizando la parte con alta densidad
Cel_Descartadas = Im1_C_Seg & ~Segmentacion_Celulas;
figure();imshow(Cel_Descartadas); title('Celulas Descartadas','Fontsize',15);
Cel_Descartadas_BW = imopen(Cel_Descartadas, B1); 
figure();imshow(Cel_Descartadas_BW); title('Celulas Descartadas','Fontsize',15);
Cel_Descartadas_Gris = Cel_Descartadas_BW.*Im1_C;
figure();imshow(Cel_Descartadas_Gris); title('Celulas Descartadas','Fontsize',15);

% Obteniendo el histograma de la imagen
[H,B] = imhist(Cel_Descartadas_Gris, 256);
H2 = H(2:end);
B2 = B(2:end);
figure();
stem(B2, H2);
title('Histograma normalizado de Fig0310(a)'); 

% Umbralizacion
Celulas2 = Cel_Descartadas_Gris > 50/255;
Celulas2_Seg = imopen(Celulas2, B1); 
% Im1_C_Seg = imfill(Im1_C_Seg,'holes');
figure();imshow(Celulas2_Seg);title('Nueva umbralizacion de las celulas recuperadas','Fontsize',15);

% Repitiendo varias veces la separacion de celulas
for i = 1:3
    % Obteniendo las propiedades de los nucleos recuperados
    [L6, n6, stats6, Tabla6] = Calcular_Propiedades(Celulas2_Seg);
    % Separando nucleos
    Condiciones_Separacion2 = [stats6.MajorAxisLength./stats6.MinorAxisLength] > 1.5 & [stats6.Area] > (area_condicion)*1.5;
    Celulas2_Seg = Separar_Nucleos(Celulas2_Seg, Condiciones_Separacion2);
    figure();imshow(Celulas2_Seg); title('Segmentacion de las celulas recuperadas','Fontsize',15);
end

Segmentacon_Final = Celulas2_Seg + Segmentacion_Celulas;
figure();imshow(Segmentacon_Final);title('Segmentacion final','Fontsize',15);


Cel_Descartadas = Im1_C_Seg & ~Segmentacion_Celulas;
[L, n] = bwlabel(Cel_Descartadas, 8);
figure(); imshow(L, []);
R5 = L == 5;
I5 = Im1_C.*R5;
figure(); imshow(I5, []);

[H,B] = imhist(I5, 256);
H2 = H(2:end);
B2 = B(2:end);
figure();
stem(B2, H2);
title('Histograma normalizado de Fig0310(a)'); 
Umbral5 = I5 > 0.4;
figure(); imshow(Umbral5, []);

g = localthresh(Im1_C, ones(9), 0, 0.9);
figure(); imshow(g, []);

[L, n] = bwlabel(g, 8);
figure(); imshow(L, []);
I5_Rellena = imfill(g,'holes');
figure(); imshow(I5_Rellena, []);

[L7, n7, stats7, Tabla7] = Calcular_Propiedades(I5_Rellena);
figure(); imshow(L7, []);

B = strel('disk', 5);
Apertura = imopen(I5_Rellena , B);
figure(); imshow(Apertura, []);

[L8, n8, stats8, Tabla8] = Calcular_Propiedades(Apertura);
figure(); imshow(L8, []);




% Imagen de las celulas enumeradas
figure();imshow(Apertura); title('Celulas enumeradas','Fontsize',15);
hold on;
for k=1:n8
    [r,c] = find(L8==k);
    rbar = mean(r);
    cbar = mean(c);
    text(cbar, rbar, num2str(k), 'Color', 'b');
end


% figure();imshow(Im1_C); title('Celulas enumeradas','Fontsize',15);
% hold on;
% for k=1:n5
%     [r,c] = find(L5==k);
%     rbar = mean(r);
%     cbar = mean(c);
%     text(cbar, rbar, num2str(k), 'Color', 'r');
% end

