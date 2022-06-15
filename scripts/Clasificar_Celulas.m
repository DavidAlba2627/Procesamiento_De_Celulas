function [Im_rgb, Tabla_Propiedades] = Clasificar_Celulas(Segmentacion_Celulas, Celulas_Duda, Segmentacion_Canal, Im_C, Im_F)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[L, n, stats, Tabla] = Calcular_Propiedades(Segmentacion_Celulas);
Num = [1:n]';
area = stats.Area;
eje_Mayor = stats.MajorAxisLength;
eje_Menor = stats.MinorAxisLength;
orientacion = stats.Orientation;
Inten_media_nucleo = cell2mat(struct2cell(regionprops(L, Im_C, 'MeanIntensity')));
Inten_media_proteina = cell2mat(struct2cell(regionprops(L, Im_F, 'MeanIntensity')));
L_Celulas_Duda = L.*Celulas_Duda;


% Calculando positivos/negativos
figure();
h_Im = histogram(L,n+1);
h_Im_val = h_Im.Values;
h_Im_val = h_Im_val(2:end);
p_o_n = zeros(length(n));
[r,c] = size(Segmentacion_Celulas);
Im_rgb = zeros(r,c,3);
porcentaje_a = zeros(n,1);
a_protenia = zeros(n,1);
for i = 1:n
     Celula = L == i;
     Celula_Duda = L_Celulas_Duda == i;
     Area_C = h_Im_val(i);
     Superposicion = Celula.*Segmentacion_Canal;
     Area_P = sum(Superposicion(:));
     
     if max(Celula_Duda(:)) == logical(0)
         if Area_P < .2*Area_C
           p_o_n(i) = 0;
           Im_rgb(:,:,3) = Im_rgb(:,:,3) + Celula;
    %      elseif .3*Area_C <= Area_P && Area_P < .6*Area_C
    %        p_o_n(i) = 1;
    %        Im_rgb(:,:,1) = Im_rgb(:,:,1) + Celula;
         else
           p_o_n(i) = 1;
           Im_rgb(:,:,2) = Im_rgb(:,:,2) + Celula;
         end
     else
         p_o_n(i) = 2;
         Im_rgb(:,:,1) = Im_rgb(:,:,1) + Celula;
     end
     
     porcentaje_a(i,1) = (Area_P*100)/Area_C;
     a_protenia(i,1) = Area_P;
end

% Im_rgb(:, :, 2) = Im_rgb(:, :, 2) - Im_rgb(:, :, 2).*Celulas_Duda;
% Im_rgb(:, :, 3) = Im_rgb(:, :, 3) - Im_rgb(:, :, 3).*Celulas_Duda;
% Im_rgb(:, :, 1) = Im_rgb(:, :, 1) + Celulas_Duda;
% Im_rgb = Im_auxG(:, :, 1) + Celulas_Duda;

imshow(Im_rgb);title('Positivos y Negativos','Fontsize',15);
if isempty(Num)
     p_o_n = [];
end
Tabla_Propiedades = table(Num, area, a_protenia, porcentaje_a, eje_Mayor, eje_Menor, orientacion, Inten_media_nucleo', Inten_media_proteina', p_o_n');
Tabla_Propiedades.Properties.VariableNames = {'Elemento','Area','Area proteina dentro de celuala','Porcentaje proteina dentro de celula', 'Eje mayor', 'Eje menor', 'Orientacion', 'Inten_media nucleo', 'Inten_media_proteina','Estado'};


a=1;
end


