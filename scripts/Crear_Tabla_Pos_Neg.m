function [CSV] = Crear_Tabla_Pos_Neg(Archivo, L, Tabla_Propiedades_A, Tabla_Propiedades_B, nombre1, nombre2, verbose)

Tabla_Propiedades_G = Tabla_Propiedades_A;
Tabla_Propiedades_IR = Tabla_Propiedades_B;
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

[r,c] = size(L);
Im_GP_IRP = zeros(r,c);
Im_GP_IRN = zeros(r,c);
Im_GN_IRP = zeros(r,c);
Im_GN_IRN = zeros(r,c);

for i=1:size(Tabla_Propiedades_G,1)
    Celula = L==i;
    Im_GP_IRP = Im_GP_IRP + Celula*GP_IRP(i); 
    Im_GP_IRN = Im_GP_IRN + Celula*GP_IRN(i); 
    Im_GN_IRP = Im_GN_IRP + Celula*GN_IRP(i); 
    Im_GN_IRN = Im_GN_IRN + Celula*GN_IRN(i);
end

CSV = Crear_Tabla_Marcajes2(Archivo, GP_IRP, GN_IRP, GN_IRN, GP_IRN, Tabla_Propiedades_A, nombre1, nombre2);

if verbose
    figure();imshow(Im_GP_IRP); title(strcat(nombre1, " Positivo - ", nombre2, " Positivo"),'Fontsize',15);
    figure();imshow(Im_GP_IRN); title(strcat(nombre1, " Positivo - ", nombre2, " Negativo"),'Fontsize',15);
    figure();imshow(Im_GN_IRP); title(strcat(nombre1, " Negativo - ", nombre2, " Positivo"),'Fontsize',15);
    figure();imshow(Im_GN_IRN); title(strcat(nombre1, " Negativo - ", nombre2, " Negativo"),'Fontsize',15);
end
end