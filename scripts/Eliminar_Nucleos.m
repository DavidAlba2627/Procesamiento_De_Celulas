function [Celulas_Depuradas] = Eliminar_Nucleos(Im_C_Seg, Condiciones)
% Esta funcion se usa para eliminar las celulas que cumplen ciertas
% condiciones

cc = bwconncomp(Im_C_Seg); 
idx = find(Condiciones); 
BW = ismember(labelmatrix(cc),idx);  
Celulas_Depuradas = Im_C_Seg -  BW;
end

