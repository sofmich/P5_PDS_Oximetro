%%======================================================================%%
%  Aplicar diversas t�cnicas del procesamiento digital de se�ales en el 
%  an�lisis y procesamiento de las se�ales proporcionadas por un Ox�metro
%  para calcular el ritmo cardiaco y el nivel de saturaci�n del oxigeno en
%  la sangre.

%%         Autor: Naum Jahaziel Nu�o Contreras
%                 Sof�a Michel Salazar Valdovinos
%%         Fecha: Mayo, 2020
%%      Profesor: Arturo Pardi�as Mir
%%    Asignatura: Procesamiento Digital de Se�ales
%%======================================================================%%

close all;
%{
    Cargar los datos de los archivos que contienen la informaci�n a
    procesar de las se�ales IR y R 
%}
load('oxi1.mat')
load('oxi2.mat')
load('oxi3.mat')