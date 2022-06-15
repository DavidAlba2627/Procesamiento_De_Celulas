function [L,n,stats,Tabla] = Calcular_Propiedades(Im_C_Seg)
% Esta funcion se usa para calcular las propiedades de las celulas 
% de una imagen usando regionprops

[L, n] = bwlabel(Im_C_Seg, 8);
stats = regionprops('table', L, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Circularity');
area = stats.Area;
eje_Mayor = stats.MajorAxisLength;
eje_Menor = stats.MinorAxisLength;
circ = stats.Circularity;
orientacion = stats.Orientation;
Num = [1:n]';
Tabla = table(Num, area, eje_Mayor, eje_Menor, orientacion, circ);
Tabla.Properties.VariableNames = {'Elemento','Area', 'Eje mayor', 'Eje menor', 'Orientation', 'Circularidad'};
end

