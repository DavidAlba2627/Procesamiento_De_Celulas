close all;
clear variables;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameetros del codigo
verbose = false;
% Este es el parametro optimizar el tamaño de nucleos entre especies de raton/ratas
area_celula = 300; %Area promedio del nucelo p/rata
%area_celula = 100; %Area promedio del nucelo p/raton

% Leyendo la imagen y separando por canal de color
cd ..; cd cell_images;
Taiep = 'Taiep2_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)'; %Nombre de los archivos se escoge entre uno u otro
SD = 'SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)';
Taiep_2 = 'TAIEP6_VentLat_TRa488-Sox2-555_DAPI_ 25x-3';
SD_2 = 'SD1_VentLat_TRa488-Sox2-555_DAPI_ 25x';
Raton = 'Raton_4m_Hp_Hoechst_SOX2488_Tra1555_TRb1647_25X_1 Zs';
Raton2 = 'Raton_4m_GD_Hoechst_Calbindin488_NeuN555_TRa1647_25X_1 Zs';
Raton3 = 'Raton_4m_Hp_Hoechst_NeuN2488_Tra1555_TRb1647_25X_1 Zs';

Archivo = Taiep;
Im1 = im2double(imread(append(Archivo,'DAPI.tif')));
cd ..; cd scripts;
%figure();imshow(Im1); title('Imagen original','Fontsize',15);
Im1_C = Im1; 
figure();imshow(Im1_C); title('Imagen Original en escale de grises','Fontsize',15);

% Realizando suavizamiento
K = fspecial('gaussian',5, 0.7);
disp('Suavizamos las celulas ...')
Im1_C_Filt = imfilter(Im1_C, K);

% Realizando el umbral global con el metodo de otsu
otsu = Im1_C_Filt < 0.2;
otsu = Im1_C_Filt.*otsu;
[h,b] = imhist(otsu, 256);
haux = h(2:end); %Eliminamos la primer componente pues esta es 0
[T,SM] = otsuthresh(haux);
Im_Seg = Im1_C_Filt > T;
Im1_C = Im_Seg.*Im1_C;


% Realizando el umbral local
disp('Aplicamos localtresh y rellenamos agujeros...')
g = localthresh(Im1_C, ones(9), 0, 0.9);
%g = localthresh(Im1_C, ones(15), 0, 0.9);
% Rellenando los agujeros
Celulas_Rellenas = imfill(g,'holes');
% Realizando la apertura
B = strel('disk', 5);
Celulas_Apertura = imopen(Celulas_Rellenas , B);

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
% Calculando propiedades
[L1, n1, stats1, Tabla1] = Calcular_Propiedades(Celulas_Inten_Correcta);

% Depurando por area pequeña
Condicion_Area = [stats1.Area] < area_celula/2;
Celulas_Area_Correcta = Eliminar_Nucleos(Celulas_Inten_Correcta, Condicion_Area);
[L2, n2, stats2, Tabla2] = Calcular_Propiedades(Celulas_Area_Correcta);
%PARA DEBUGGING
if verbose
    figure();imshow(Celulas_Area_Correcta); title('Celulas enumeradas depuaradas en area','Fontsize',15);
    hold on;
    for k=1:n2
        [r,c] = find(L2==k);
        rbar = mean(r);
        cbar = mean(c);
        text(cbar, rbar, num2str(k), 'Color', 'b');
    end
end
% Dividiendo las celulas
Celulas_Divididas = Celulas_Area_Correcta;


[L3, n3, stats3, Tabla3] = Calcular_Propiedades(Celulas_Divididas);
% Separando nucleos
Condiciones_Separacion2 = [stats3.MajorAxisLength./stats3.MinorAxisLength] > 1.4 & [stats3.Area] > area_celula*1.5 & [stats3.Area] < area_celula*4.0 & [stats3.Circularity] > 0.4;
%Agregar parametro VERBOSE para mostrar imagenes o no
Celulas_Divididas = Separar_Nucleos(Celulas_Divididas, Condiciones_Separacion2, verbose);

% Separando las celulas en duda
[L4, n4, stats4, Tabla4] = Calcular_Propiedades(Celulas_Divididas);
Condiciones_Marcado_Rojo = [stats4.Area] > area_celula*2 & [stats4.Circularity] < 0.75;
Celulas_Correctas_Seg = Eliminar_Nucleos(Celulas_Divididas, Condiciones_Marcado_Rojo);
Celulas_Duda = Celulas_Divididas - Celulas_Correctas_Seg;
figure();imshow(Celulas_Duda.*Im1_C);title('Celulas en duda','Fontsize',15);

hold off;
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
%%%%%%%%%%%%%%%%%%%%%%%
k_aux = fspecial('gaussian',100, 50);
Im_partes = imfilter(Celulas_Divididas, k_aux, 'symmetric');
B = strel('disk', 10);
[h,b] = imhist(Im_partes, 256);
[T,SM] = otsuthresh(h(2:end));
Im_Seg = im2bw(Im_partes,T);


[Lseg, nseg, statsseg, Tablaseg] = Calcular_Propiedades(Im_Seg);
i = find(max(statsseg.Area)==statsseg.Area);
se = strel('line', statsseg.MajorAxisLength(i)*0.1, statsseg.Orientation(i)); 
Im_dilatada= imdilate(Im_Seg, se);
cerradura = strel('disk',20);
Im_dilatada = imclose(Im_dilatada,cerradura);
[Lseg2, nseg2, statsseg2, Tablaseg2] = Calcular_Propiedades(Im_dilatada);
Region1 = Lseg2 == find(statsseg2.Area==max(statsseg2.Area));
Region1Grises = Region1.*Im1_C;

Region2 = Im_partes.*~Region1;

Region_sobrante = ~Region1;
[Lregiones, nregiones] = bwlabel(Region_sobrante,8);
Region2 = Lregiones == 1;
Region2Grises = Region2.*Im1_C;
Region3 = Lregiones == 2;
Region3Grises = Region3.*Im1_C;
Region_Total = Region1 + 2*Region2 + 3*Region3;

Celulas_completas = Region1 | Celulas_Divididas;
[Lcomp, ncomp, statscomp, Tablacomp] = Calcular_Propiedades(Celulas_completas);
Region1_Corregida = Lcomp == find(statscomp.Area==max(statscomp.Area));
Region1_Corregida = imfill(Region1_Corregida,'holes');
%Esta para VERBOSE
if verbose
    figure();imshow(Region1_Corregida); title('Región central','Fontsize',15);
end


Region_sobrante = ~Region1_Corregida;
[Lregiones, nregiones] = bwlabel(Region_sobrante,8);
Region2_Corregida = Lregiones == 1;
Region2Grises_Corregida = Region2_Corregida.*Im1_C;
Region3_Corregida = Lregiones == 2;
Region3Grises_Corregida = Region3_Corregida.*Im1_C;
%Para VERBOSE estas 3
if verbose
    figure();imshow(~(Region2_Corregida+Region3_Corregida).*Im1_C,[]); title('Region 1','Fontsize',15);
    figure();imshow(Region2Grises_Corregida,[]); title('Region 2','Fontsize',15);
    figure();imshow(Region3Grises_Corregida,[]); title('Region 3','Fontsize',15);
end


cd ..; cd cell_images;
Canales = open(append(Archivo,'Canales.mat')).Indices_Nombre;
Canales_Sin_DAPI = open(append(Archivo,'Canales.mat')).Indices_Nombre_Sin_DAPI;
Escala = open(append(Archivo,'Canales.mat')).Escala;
cd ..; cd scripts;

Region_Total_Corregida = Region1_Corregida + 2*Region2_Corregida + 3*Region3_Corregida;
[stats_regiones, nucleos_p_area, area_region_i, nucleos_p_region, indices_regiones] = Calcular_Propiedades_Regiones(Celulas_Divididas, Region_Total_Corregida, Escala, L4, n4);
Tabla_Regiones = table(['Region 1';'Region 2'; 'Region 3'],nucleos_p_region, area_region_i,nucleos_p_area);
Tabla_Regiones.Properties.VariableNames = {'Regiones','Nucleos por region', 'Area (micras cuadradas)', 'Densidad (nucleo por micra cuadrada)'};

Tabla4.('Region de pertenencia') = indices_regiones;
figure();imshow(Region_Total_Corregida,[]); title('Regiones','Fontsize',15);
hold on;
for k=1:3
    [R,C] = find(Region_Total_Corregida==k);
    rbar = mean(R);
    cbar = mean(C);
    text(cbar, rbar, num2str(k), 'Color', 'r', 'Fontsize',30);
end
%%%%%%%%%%%%%%%%%%%%%%%


% Segmentacion de canales de color
% Aplicando el métod de Otsu
%Aplicar VERBOSE a las funciones
Imagenes_Canales_Seg = zeros(r,c,size(Canales_Sin_DAPI, 1));
Imagenes_Canales_Gris = zeros(r,c,size(Canales_Sin_DAPI,1));
Nombre_Canales = string(zeros(3,1));
% for i=1:size(Canales_Sin_DAPI,1)
%     [Imagenes_Canales_Seg(:,:,i),Imagenes_Canales_Seg(:,:,i)] = Seg_canales(append(Archivo,Canales_Sin_DAPI,'.tif'),2, verbose);
% end
a = 4;
if size(Canales_Sin_DAPI,1) > 3
    error('No se puede tener mas de 4 canales');
elseif size(Canales_Sin_DAPI,1) == 3
    [Im_Canal1_seg, Im_Canal1_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(1),'.tif'),1, verbose);
    [Im_Canal2_seg, Im_Canal2_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(2),'.tif'),1, verbose);
    [Im_Canal3_seg, Im_Canal3_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(3),'.tif'),1, verbose);
    nombres = Canales_Sin_DAPI;
elseif size(Canales_Sin_DAPI,1) == 2
    [Im_Canal1_seg, Im_Canal1_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(1),'.tif'),1, verbose);
    [Im_Canal2_seg, Im_Canal2_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(2),'.tif'),1, verbose);
    [Im_Canal3_seg, Im_Canal3_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(2),'.tif'),1, verbose);
    nombres = [Canales_Sin_DAPI(1); Canales_Sin_DAPI(2); Canales_Sin_DAPI(2)];
elseif size(Canales_Sin_DAPI,1) == 1
    [Im_Canal1_seg, Im_Canal1_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(1),'.tif'),1, verbose);
    [Im_Canal2_seg, Im_Canal2_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(1),'.tif'),1, verbose);
    [Im_Canal3_seg, Im_Canal3_gris] = Seg_canales(append(Archivo,Canales_Sin_DAPI(1),'.tif'),1, verbose);
    nombres = [Canales_Sin_DAPI(1); Canales_Sin_DAPI(1); Canales_Sin_DAPI(1)];
end

% Agregar titulo a la imagen de la funcion
[Im_Class_1, Tabla_Propiedades_Marcaje_1] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, Im_Canal1_seg,  Im1_C, Im_Canal1_gris);
title(strcat("Clasificacion sobre el marcaje ",nombres(1)),'Fontsize',15);
[Im_Class_2, Tabla_Propiedades_Marcaje_2] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, Im_Canal2_seg,  Im1_C, Im_Canal2_gris);
title(strcat("Clasificacion sobre el marcaje ",nombres(2)),'Fontsize',15);
[Im_Class_3, Tabla_Propiedades_Marcaje_3] = Clasificar_Celulas(Celulas_Divididas, Celulas_Duda, Im_Canal3_seg,  Im1_C, Im_Canal3_gris);
title(strcat("Clasificacion sobre el marcaje ",nombres(3)),'Fontsize',15);

% % Crear tablas para c/u de las regiones
% % Region 1
% [Im_Class_1_region_1, Tabla_Propiedades_Marcaje_1_region_1] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==1, Celulas_Duda.*Region_Total_Corregida==1, Im_Canal1_seg.*Region_Total_Corregida==1,  Im1_C.*Region_Total_Corregida==1, Im_Canal1_gris.*Region_Total_Corregida==1);
% title(strcat("Clasificacion sobre el marcaje ",nombres(1)),'Fontsize',15);
% [Im_Class_2_region_1, Tabla_Propiedades_Marcaje_2_region_1] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==1, Celulas_Duda.*Region_Total_Corregida==1, Im_Canal2_seg.*Region_Total_Corregida==1,  Im1_C.*Region_Total_Corregida==1, Im_Canal2_gris.*Region_Total_Corregida==1);
% title(strcat("Clasificacion sobre el marcaje ",nombres(2)),'Fontsize',15);
% [Im_Class_3_region_1, Tabla_Propiedades_Marcaje_3_region_1] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==1, Celulas_Duda.*Region_Total_Corregida==1, Im_Canal3_seg.*Region_Total_Corregida==1,  Im1_C.*Region_Total_Corregida==1, Im_Canal3_gris.*Region_Total_Corregida==1);
% title(strcat("Clasificacion sobre el marcaje ",nombres(3)),'Fontsize',15);
% if size(Tabla_Propiedades_Marcaje_1_region_1,1)~=0 & size(Tabla_Propiedades_Marcaje_2_region_1,1)~=0 & size(Tabla_Propiedades_Marcaje_3_region_1,1)~=0
%     Conteo_region_1 = Crear_Tabla_Marcajes(Archivo, Tabla_Propiedades_Marcaje_1_region_1,Tabla_Propiedades_Marcaje_2_region_1, Tabla_Propiedades_Marcaje_3_region_1, nombres);
% else
%     fprintf('En la region 1 no hay nucleos para realizar el conteo');
% end
% 
% %Region 2
% n=2;
% [Im_Class_1_region_2, Tabla_Propiedades_Marcaje_1_region_2] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==2, Celulas_Duda.*Region_Total_Corregida==2, Im_Canal1_seg.*Region_Total_Corregida==2,  Im1_C.*Region_Total_Corregida==2, Im_Canal1_gris.*Region_Total_Corregida==2);
% title(strcat("Clasificacion sobre el marcaje ",nombres(1)),'Fontsize',15);
% [Im_Class_2_region_2, Tabla_Propiedades_Marcaje_2_region_2] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==n, Celulas_Duda.*Region_Total_Corregida==n, Im_Canal2_seg.*Region_Total_Corregida==n,  Im1_C.*Region_Total_Corregida==n, Im_Canal2_gris.*Region_Total_Corregida==n);
% title(strcat("Clasificacion sobre el marcaje ",nombres(2)),'Fontsize',15);
% [Im_Class_3_region_2, Tabla_Propiedades_Marcaje_3_region_2] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==n, Celulas_Duda.*Region_Total_Corregida==n, Im_Canal3_seg.*Region_Total_Corregida==n,  Im1_C.*Region_Total_Corregida==n, Im_Canal3_gris.*Region_Total_Corregida==n);
% title(strcat("Clasificacion sobre el marcaje ",nombres(3)),'Fontsize',15);
% if size(Tabla_Propiedades_Marcaje_1_region_2,1)~=0 & size(Tabla_Propiedades_Marcaje_2_region_2,1)~=0 & size(Tabla_Propiedades_Marcaje_3_region_2,1)~=0
%     Conteo_region_2 = Crear_Tabla_Marcajes(Archivo, Tabla_Propiedades_Marcaje_1_region_2,Tabla_Propiedades_Marcaje_2_region_2, Tabla_Propiedades_Marcaje_3_region_2, nombres);
% else
%     fprintf('En la region 2 no hay nucleos para realizar el conteo');
% end
% 
% %Region 3
% n=3;
% [Im_Class_1_region_3, Tabla_Propiedades_Marcaje_1_region_3] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==n, Celulas_Duda.*Region_Total_Corregida==n, Im_Canal1_seg.*Region_Total_Corregida==n,  Im1_C.*Region_Total_Corregida==n, Im_Canal1_gris.*Region_Total_Corregida==n);
% title(strcat("Clasificacion sobre el marcaje ",nombres(1)),'Fontsize',15);
% [Im_Class_2_region_3, Tabla_Propiedades_Marcaje_2_region_3] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==n, Celulas_Duda.*Region_Total_Corregida==n, Im_Canal2_seg.*Region_Total_Corregida==n,  Im1_C.*Region_Total_Corregida==n, Im_Canal2_gris.*Region_Total_Corregida==n);
% title(strcat("Clasificacion sobre el marcaje ",nombres(2)),'Fontsize',15);
% [Im_Class_3_region_3, Tabla_Propiedades_Marcaje_3_region_3] = Clasificar_Celulas(Celulas_Divididas.*Region_Total_Corregida==n, Celulas_Duda.*Region_Total_Corregida==n, Im_Canal3_seg.*Region_Total_Corregida==n,  Im1_C.*Region_Total_Corregida==n, Im_Canal3_gris.*Region_Total_Corregida==n);
% title(strcat("Clasificacion sobre el marcaje ",nombres(3)),'Fontsize',15);
% if size(Tabla_Propiedades_Marcaje_1_region_3,1)~=0 & size(Tabla_Propiedades_Marcaje_2_region_3,1)~=0 & size(Tabla_Propiedades_Marcaje_3_region_3,1)~=0
%     Conteo_region_3 = Crear_Tabla_Marcajes(Archivo, Tabla_Propiedades_Marcaje_1_region_3,Tabla_Propiedades_Marcaje_2_region_3, Tabla_Propiedades_Marcaje_3_region_3, nombres);
% else
%     fprintf('En la region 3 no hay nucleos para realizar el conteo');
% end


% Tabla_Propiedades_G = Tabla_Propiedades_1;
% Tabla_Propiedades_IR = Tabla_Propiedades_2;
% ID_G = Tabla_Propiedades_G(:, 'Estado');
% ID_G = ID_G{:,:};
% ID_IR = Tabla_Propiedades_IR(:, 'Estado');
% ID_IR = ID_IR{:,:};
% 
% % Etiquetas de los canales verde e infrarrojo
% G_Pos = ID_G == 1;
% G_Neg = ID_G == 0;
% IR_Pos = ID_IR == 1;
% IR_Neg = ID_IR == 0;
% % Ambos positivos
% GP_IRP = G_Pos & IR_Pos;
% % Verde positivo e infrarrojo negativo
% GP_IRN = G_Pos & IR_Neg;
% % Verde negativo e infrarrojo positivo
% GN_IRP = G_Neg & IR_Pos;
% % Ambos negativos
% GN_IRN = G_Neg & IR_Neg;
% 
% Im_GP_IRP = zeros(r,c);
% Im_GP_IRN = zeros(r,c);
% Im_GN_IRP = zeros(r,c);
% Im_GN_IRN = zeros(r,c);
% 
% for i=1:n4
%     Celula = L4==i;
%     Im_GP_IRP = Im_GP_IRP + Celula*GP_IRP(i); 
%     Im_GP_IRN = Im_GP_IRN + Celula*GP_IRN(i); 
%     Im_GN_IRP = Im_GN_IRP + Celula*GN_IRP(i); 
%     Im_GN_IRN = Im_GN_IRN + Celula*GN_IRN(i);
% end
% 
% 
% figure();imshow(Im_GP_IRP); title('Verde positivo e infrarrojo positivo','Fontsize',15);
% figure();imshow(Im_GP_IRN); title('Verde positivo e infrarrojo negativo','Fontsize',15);
% figure();imshow(Im_GN_IRP); title('Verde negativo e infrarrojo positivo','Fontsize',15);
% figure();imshow(Im_GN_IRN); title('Verde negativo e infrarrojo negativo','Fontsize',15);
% figure();imshow(Im_GP_IRP + Im_GP_IRN + Im_GN_IRP + Im_GN_IRN); title('Verde negativo e infrarrojo negativo','Fontsize',15);

%Agregando columna de region de pertenencia
Tabla_Propiedades_Marcaje_1.('Region de pertenencia') = indices_regiones;
Tabla_Propiedades_Marcaje_2.('Region de pertenencia') = indices_regiones;
Tabla_Propiedades_Marcaje_3.('Region de pertenencia') = indices_regiones;


%Creando un .csv para los canales
Conteo_general = Crear_Tabla_Marcajes(Archivo, Tabla_Propiedades_Marcaje_1,Tabla_Propiedades_Marcaje_2, Tabla_Propiedades_Marcaje_3, nombres);
% CSV2 = Crear_Tabla_Marcajes2(Archivo, GP_IRP, GN_IRP, GN_IRN, GP_IRN, Tabla_Propiedades_3);
Expresion_canal1_canal2 = Crear_Tabla_Pos_Neg(Archivo, L4, Tabla_Propiedades_Marcaje_1, Tabla_Propiedades_Marcaje_2, nombres(1), nombres(2), verbose);
Expresion_canal1_canal3 = Crear_Tabla_Pos_Neg(Archivo, L4, Tabla_Propiedades_Marcaje_1, Tabla_Propiedades_Marcaje_3, nombres(1), nombres(3), verbose);
Expresion_canal2_canal3 = Crear_Tabla_Pos_Neg(Archivo, L4, Tabla_Propiedades_Marcaje_2, Tabla_Propiedades_Marcaje_3, nombres(2), nombres(3), verbose);