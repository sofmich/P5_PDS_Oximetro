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
    procesar de las se�ales IR y R, as� como la respectiva frecuencia de 
    muestreo
%}

%% Lectura de los datos
oxi1 = load('oxi1.mat','-mat');
oxi2 = load('oxi2.mat','-mat');
oxi3 = load('oxi3.mat','-mat');
% igualamos el valor de oxi actual al leido
oxi_actual = oxi3;
%eje para graficar
n=1:size(oxi_actual.x_ir);
%transpuesta
oxi_actual_ir=oxi_actual.x_ir.';
oxi_actual_red=oxi_actual.x_red.';
%frecuencia de muestreo
Fs=oxi_actual.fs;
%eje de tiempo
dt=1/Fs;
t=0:dt:length(oxi_actual_ir)/(Fs);
t=t(2:length(oxi_actual_ir)+1);

%% graficamos lo leido

figure;plot(t,oxi_actual_ir );
title('Lectura infrarojo oxi1');
figure;plot(t,oxi_actual_red );
title('Lectura red oxi1');
%% Vectores x para graficar frecuencias

oxi_red=oxi_actual_red;
oxi_ir=oxi_actual_ir;

fourier_lenght=length(oxi_ir);
%rango de frecuencia rads
w=0:2*pi/(fourier_lenght-1):2*pi;
fre=(w*Fs)/(2*pi);

%% Dise�o del filtro passa bandas con ventana triangular para eliminar rizo
%60 ppm
wp = 0.5/ (Fs/2);
%230 ppm
wr= 3.8/(Fs/2);
% Filtro pasa bandas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bandpass = fir1(100,[wp wr], 'bandpass', triang(101));
freqz(bandpass); title('filtro pasa bandas');
%% Eliminar componente directa de los datos num�ricos de R, IR
ir = oxi_ir - mean(oxi_ir);
red = oxi_red - mean(oxi_red);

%% Apllicar el filtro a las nuevas se�ales sin componente directa
p_banda_red=filter(bandpass,1,red);
p_banda_ir=filter(bandpass,1,ir);

%% Sacamos promedio de max y mins de la señal 
%% C�lculo de saturaci�n de ox�geno en la sangre
vals_max_oxi1_red = abs(findpeaks(p_banda_red));
average_max_oxi1_red = mean(vals_max_oxi1_red);

% Findpeaks encuentra los m�ximos relativos, para su uso es necesario 
% hacer un promedio 
vals_max_oxi1_ir = abs(findpeaks(p_banda_ir));
average_max_oxi1_ir = mean(vals_max_oxi1_ir);

%% Encontrar los valores m�nimos de las se�ales, R e IR
vals_min_oxi1_red = abs(findpeaks(-p_banda_red));
average_min_oxi1_red = mean(vals_min_oxi1_red);

vals_min_oxi1_ir = abs(findpeaks(-p_banda_ir));
average_min_oxi1_ir = mean(vals_min_oxi1_ir);

%% Componente en AC para poder realizar las operaciones matem�ticas
AC_oxi1_red = abs(average_max_oxi1_red - average_min_oxi1_red);
AC_oxi1_ir = abs(average_max_oxi1_ir- average_min_oxi1_ir);

%% Graficar la transformada de Fourier de las se�ales sin componente 
%% directa y filtradas 

FFT_IR = fft(p_banda_ir);
%graficas de las TFF
figure; plot(fre,abs(FFT_IR));
title('Espectro Frecuencias infrarojo oxi1');

%%FFT for RED signal
FFT_RED = fft(p_banda_red);
%graficas de las TFF
figure;plot(fre,abs(FFT_RED));
title('Espectro Frecuencias rojo oxi1');
%%  Encontrar el maximo local en el espectro rojo para calcular los cambios
%%  entre picos de la se�al y encontrar los latiduos por minuto

% M�ximo de la se�al de luz Red, el espectro permite identificar la freq.
% fundamental que indica el elemento con una mayor magnitud 
[max_MAGR, max_RED]= max(abs(FFT_RED));
max_RED = fre(max_RED);
% M�ximo de la se�al de luz Infrarrojo, el espectro permite identificar 
% la freq. fundamental que indica el elemento con una mayor magnitud 
[max_MAGIR, max_IR] = max(abs(FFT_IR));% 
max_IR = fre(max_IR);

%% Calcular los latidos por minuto realizando promedio entre los m�ximos
%% y considerando que 1 minuto = 60 pulsaciones (1 cada segundo)
HEART_RATE_oxi1 = mean([max_IR, max_RED])*60;

%% Calcular la concentraci�n de ox�geno en la sangre aplicando un logaritmo
%% a la se�al, esto se puede encontrar en documentaci�n de diferentes f
%% fuentes donde se aplica una formula: R = log(ACred) / log(ACir) una
%% relaci�n entre las componentes directas y de ac de los valores del oxim.
% Calculo de Sp02
R_RED1 =log(AC_oxi1_red);

R_IR1=log(AC_oxi1_ir);

R1=R_RED1/R_IR1;
%% El siguiente paso es realizar el calculo con f�rmula 
%% C�lculado a trav�s de la formula de la comunidad de matlab
%% SpO2 = 100 - 25 * R 
SpO2_oxi1 = 110 - 25*R1;



