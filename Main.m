clear variables; close all; clc;
% Paramametros que el usuario debe ingresar

%Colocar el nombre del archivo (sin extension)
%Se asume que existe el archivo .czi y está colocado 
%en la carpeta Color_channels/Images
file_name = 'SD1_6m_PLV_BrdU488_GFAP555_SOX2647_25X (1)';
%El verbose en estado True es para mostrar imagenes de pasos intermedios
verbose = false;


% Funciones que realizan el conteo
%Las siguientes lineas realizan una separacion en canales de la imagen
%original con extension .czi, generando archivos en formato .tif, estos se
%guardan en la carpeta Cell_images
cd Color_channels;
FuncionSepararCanales(file_name);

%Las siguientes lineas realizan el proceso de analis de las imagenes
%Los resultados se guardan en la carpeta Tablas_de_estadisticas
cd ..; cd scripts;
ContarCelulas(file_name,verbose);