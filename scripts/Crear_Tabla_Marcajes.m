function [CSV] = Crear_Tabla_Marcajes(Archivo, Tabla_Propiedades_G,Tabla_Propiedades_R, Tabla_Propiedades_IR, nombres)
Etiquetas_G = Tabla_Propiedades_G{:,10};
[C,ia,ic_G] = unique(Etiquetas_G);
Conteo_Etiquetas_G = accumarray(ic_G,1);

Etiquetas_R = Tabla_Propiedades_R{:,10};
[C,ia,ic_R] = unique(Etiquetas_R);
Conteo_Etiquetas_R = accumarray(ic_R,1);

Etiquetas_IR = Tabla_Propiedades_IR{:,10};
[C,ia,ic_IR] = unique(Etiquetas_IR);
Conteo_Etiquetas_IR = accumarray(ic_IR,1);

Conteo = vertcat(Conteo_Etiquetas_G, Conteo_Etiquetas_R, Conteo_Etiquetas_IR);
Marcaje = [[nombres(1), nombres(1), nombres(1)]'; [nombres(2), nombres(2), nombres(2)]'; [nombres(3), nombres(3), nombres(3)]'];
Etiqueta = [["Negativo", "Positivo", "Duda"]';["Negativo", "Positivo", "Duda"]';["Negativo", "Positivo", "Duda"]'];

CSV = table(Marcaje, Etiqueta, Conteo);
CSV.Properties.VariableNames = {'Marcaje','Etiqueta', 'Conteo'};

cd ..; cd tables_graphics;
writetable(CSV, strcat('tabla1_', Archivo,'.csv'));
cd ..; cd scripts;
end

