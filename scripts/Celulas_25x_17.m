close all;
clear variables;
clc;

% Leyendo la imagen y separando por canal de color
cd ..; cd cell_images;
Im1 = im2double(imread('25x-NeuN-Hoechst-izq-17 ZEN.tif'));
cd ..; cd scripts;
figure();imshow(Im1); title(' Imagen original','Fontsize',15);
Im1_C = Im1(:,:,3); 
Im1_F = Im1(:,:,2);
figure();imshow(Im1_C); title(' Imagen de las celulas','Fontsize',15);
figure();imshow(Im1_F); title(' Imagen de la proteina','Fontsize',15);

% Segmentando los nucleos
Celulas = Im1_C > 15/255;
B1 = strel('disk', 5);  
Im1_C_Seg = imopen(Celulas, B1); 
Im1_C_Seg = imfill(Im1_C_Seg,'holes');
% figure();imshow(Im1_C_Seg);title('1. Segmentacion de las celulas','Fontsize',15);

% Obteniendo las propiedades de los nucleos
[L1, n1, stats1, Tabla1] = Calcular_Propiedades(Im1_C_Seg);
% Eliminando nucleos que no nos importan (Grandes y no circulares)
area_condicion = 270;
Condiciones = [stats1.Area] < area_condicion/3 | [stats1.Area] > area_condicion*3 | [stats1.Circularity] < 0.3;
Celulas_Depuradas = Eliminar_Nucleos(Im1_C_Seg, Condiciones);
figure();imshow(Celulas_Depuradas); title('Celulas buenas','Fontsize',15);

% Obteniendo las propiedades de los nucleos depurados
[L2, n2, stats2, Tabla2] = Calcular_Propiedades(Celulas_Depuradas);
% Separando nucleos
Condiciones_Separacion = [stats2.MajorAxisLength./stats2.MinorAxisLength] > 1.5 & [stats2.Area] > (area_condicion-200)*1.5;
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
figure(); imshow(Celulas_Crecidas); title('Imagen de crecimiento de regiones','Fontsize',15);
% Eliminando celulas que crecieron mucho
[L4, n4, stats4, Tabla4] = Calcular_Propiedades(Celulas_Crecidas);
Condiciones_Crecidas = [stats4.Area] < area_condicion/3 | [stats4.Area] > area_condicion*2 | [stats4.Circularity] < 0.3;
Celulas_No_Crecidas = Eliminar_Nucleos(Celulas_Crecidas, Condiciones_Crecidas);
Segmentacion_Celulas = Celulas_No_Ovaladas + Celulas_No_Crecidas; 
figure(); imshow(Segmentacion_Celulas); title('Segmemtacion de celulas','Fontsize',15);

% ------ Proteina (antes farmaco) ---------
Im1_F_Seg = Im1_F > 0.14;
figure();imshow(Im1_F_Seg, []); title('Segmentacion de la proteina','Fontsize',15);

% Tabla de propiedades de las celulas
[L5, n5, stats5, Tabla5] = Calcular_Propiedades(Segmentacion_Celulas);
Num5 = [1:n5]';
area5 = stats5.Area;
eje_Mayor5 = stats5.MajorAxisLength;
eje_Menor5 = stats5.MinorAxisLength;
orientacion5 = stats5.Orientation;
Inten_media_nucleo = cell2mat(struct2cell(regionprops(L5, Im1_C, 'MeanIntensity')));
Inten_media_proteina = cell2mat(struct2cell(regionprops(L5, Im1_F, 'MeanIntensity')));

% Calculando positivos/negativos
figure();
h_Im = histogram(L5,n5+1);
h_Im_val = h_Im.Values;
h_Im_val = h_Im_val(2:end);
p_o_n = zeros(length(n5));
[r,c] = size(L4);
Im_rgb = zeros(r,c,3);
porcentaje_a = zeros(n5,1);
a_protenia = zeros(n5,1);
for i = 1:n5
     Celula = L5 == i;
     Area_C = h_Im_val(i);
     Superposicion = Celula.*Im1_F_Seg;
     Area_P = sum(Superposicion(:));
     
     if Area_P < .3*Area_C
       p_o_n(i) = 0;
       Im_rgb(:,:,3) = Im_rgb(:,:,3) + Celula;
     else
       p_o_n(i) = 1;
       Im_rgb(:,:,2) = Im_rgb(:,:,2) + Celula;
     end
     porcentaje_a(i,1) = (Area_P*100)/Area_C;
     a_protenia(i,1) = Area_P;
end

figure();imshow(Im_rgb);title('Positivos y Negativos','Fontsize',15);
Tabla_Propiedades = table(Num5, area5, a_protenia, porcentaje_a, eje_Mayor5, eje_Menor5, orientacion5, Inten_media_nucleo', Inten_media_proteina', p_o_n');
Tabla_Propiedades.Properties.VariableNames = {'Elemento','Area','Area proteina dentro de la celuala','Porcentaje', 'Eje mayor', 'Eje menor', 'Orientacion', 'Inten_media nucleo', 'Inten_media_proteina','Estado'};


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

% Imprimiendo un subplot
figure();subplot(1,2,1);imshow(Im1);title('Imagen original','Fontsize',15);
subplot(1,2,2);imshow(Im_rgb);title("Positivos y Negativos",'Fontsize',15);