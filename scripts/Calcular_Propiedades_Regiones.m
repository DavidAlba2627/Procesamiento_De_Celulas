function [stats_regiones, nucleos_p_area, area_region_i, nucleos_p_region, indices_regiones] = Calcular_Propiedades_Regiones(img_original,img_regiones, escala, L4, n4)

%La funcion encuentra la cantidad de regiones presentes en la imagen de
%regiones. Calcula las propiedades de c/u de las celulas de c/u de las
%regiones.
%Esto lo guarda en un struct de tamanio de la cant. de regiones, donde cada
%espacio es la tabla correspondinte. Y cada instancia de las tablas es un
%nucleo.

n_regiones = max(img_regiones(:));
stats_regiones = cell(n_regiones, 1);

%Loop para c/u de los nucleos de c/u de las regiones
nucleos_p_area = zeros(n_regiones, 1); %esto para proporcion area/nucleos
nucleos_p_region = zeros(n_regiones, 1);
area_region_i = zeros(n_regiones, 1);
indices_regiones = zeros(n4, 1);
for i = 1:n_regiones
    img_aux = immultiply(img_original,(img_regiones==i));
    %figure(); imshow(img_aux);
    [L, n] = bwlabel(img_aux, 8);
    stats_regiones{i} = regionprops('table', L, img_original, 'Area', 'MeanIntensity', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Circularity');
    %Seccion para la proporcion area/nucleos
    area_region_i(i) = cell2mat(struct2cell(regionprops(img_regiones == i, 'Area')));
    area_region_i(i) = escala(1) * escala(2) * cell2mat(struct2cell(regionprops(img_regiones == i, 'Area')));
    nucleos_region_i = size(stats_regiones{i}, 1);
    nucleos_p_area(i) = nucleos_region_i /   area_region_i(i);
    nucleos_p_region(i) = nucleos_region_i;
end


for i = 1:n4
   aux = L4==i;
   for j = 1:n_regiones
       if max(max(aux.*(img_regiones==j))) == 1
           indices_regiones(i) = j;
       end
   end
end
end