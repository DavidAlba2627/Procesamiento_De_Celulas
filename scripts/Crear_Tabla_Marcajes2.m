function [CSV] = Crear_Tabla_Marcajes2(Archivo, GP_IRP, GN_IRP, GN_IRN, GP_IRN, Tabla_Propiedades_G, nombre1, nombre2)

Etiquetas_G = Tabla_Propiedades_G{:,10};
[C,ia,ic_G] = unique(Etiquetas_G);
Conteo_Etiquetas_G = accumarray(ic_G,1);


Marcaje =  ["Total", strcat(nombre1,"+/", nombre2,"+"), strcat(nombre1,"-/",nombre2,"+"),strcat(nombre1,"-/",nombre2,"-"),strcat(nombre1,"+/",nombre2,"-"), "Warning"]';
%Marcaje =  ["Total", "Verde+/IR+", "Verde-/IR+", "Verde-/IR-", "Verde+/IR-", "Warning"]';
Conteo = [size(GP_IRN,1), sum(GP_IRP), sum(GN_IRP), sum(GN_IRN), sum(GP_IRN), Conteo_Etiquetas_G(end)]';
CSV = table(Marcaje, Conteo);
CSV.Properties.VariableNames = {'Marcaje', 'Conteo'};

cd ..; cd tables_graphics;
writetable(CSV, strcat('tabla2_', Archivo,'.csv'));
cd ..; cd scripts;
end
