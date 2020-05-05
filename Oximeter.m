%%======================================================================%%
%  Aplicar diversas técnicas del procesamiento digital de señales en el 
%  análisis y procesamiento de las señales proporcionadas por un Oxímetro
%  para calcular el ritmo cardiaco y el nivel de saturación del oxigeno en
%  la sangre.

%%         Autor: Naum Jahaziel Nuño Contreras
%                 Sofía Michel Salazar Valdovinos
%%         Fecha: Mayo, 2020
%%      Profesor: Arturo Pardiñas Mir
%%    Asignatura: Procesamiento Digital de Señales
%%======================================================================%%

close all;
%{
    Cargar los datos de los archivos que contienen la información a
    procesar de las señales IR y R 
%}
load('oxi1.mat')
load('oxi2.mat')
load('oxi3.mat')