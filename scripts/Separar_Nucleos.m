function [Celulas_Separadas] = Separar_Nucleos(Celulas_Depuradas, Condiciones, verbose)
% Esta funcion se usa para separar los nucleos donde parece que se unieron
% 2 celulas diferentes

% Seleccionando los nucleos unidos
cc = bwconncomp(Celulas_Depuradas); 
idx = find(Condiciones); 
BW = ismember(labelmatrix(cc),idx);
if verbose
    figure();imshow(BW); title('Posibles celulas unidas','Fontsize',15);
end
Celulas_Sin_Pegar = Celulas_Depuradas -  BW;

% Dividiendo las celulas
[L, n] = bwlabel(BW, 8);
stats = regionprops('table', L, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
eje_Menor = stats.MinorAxisLength;
orientacion = stats.Orientation;
Area = stats.Area;
[r,c] = size(Celulas_Depuradas);
Celulas_Divididas = zeros(r,c);
B = strel('disk', 3);  

for i = 1:n
    Nucleo = L == i;
    se2 = strel('line', eje_Menor(i), orientacion(i)+90); 
    n_aux = 1;
    k = 0;
    while n_aux == 1
        se1 = strel('line', eje_Menor(i)*(0.25+0.02*k), orientacion(i)+90);
        k = k + 1;
        Im_Segmentacion = imerode(Nucleo, se1);
        %figure();imshow(Im_Segmentacion); title('Ovaladas','Fontsize',15);
        Im_Segmentacion = imopen(Im_Segmentacion, B); 
        %figure();imshow(Im_Segmentacion); title('Ovaladas','Fontsize',15);
        Im_Segmentacion = imdilate(Im_Segmentacion, se2);
        %figure();imshow(Im_Segmentacion); title('Celula','Fontsize',15);
        Celula_Dividida = Nucleo.*Im_Segmentacion;
        [L_aux, n_aux, stats_aux, Tabla_aux] = Calcular_Propiedades(Celula_Dividida);
        Area_aux = stats_aux.Area;
        if n_aux == 0
            Celula_Dividida = Nucleo;
        end
    end
    
    if sum(Area_aux) < 0.6*Area(i)
        Celula_Dividida = Nucleo;    
    end
    
    Celulas_Divididas = Celulas_Divididas + Celula_Dividida;
end

Celulas_Separadas = Celulas_Sin_Pegar + Celulas_Divididas;
if verbose
    figure();imshow(Celulas_Divididas);title('Division de los posibles nucleos unidos','Fontsize',15);
end
end

