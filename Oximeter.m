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
%% Alternar entre los distintos conjutnos de datos
oxi_actual = oxi3;
% Eje horizontal para graficar
n=1:size(oxi_actual.x_ir);
% Transpuesta para poder manipular los vectores con facilidad
oxi_actual_ir=oxi_actual.x_ir.';
oxi_actual_red=oxi_actual.x_red.';
% Frecuencia de muestreo de los datos para poder procesar
Fs=oxi_actual.fs;
% Eje de tiempo
dt=1/Fs;
t=0:dt:length(oxi_actual_ir)/(Fs);
t=t(2:length(oxi_actual_ir)+1);

%% Grafica de los datos le�dos
figure;plot(t,oxi_actual_ir );
title('Lectura infrarojo (IR)');
ylabel('Amplitud');
xlabel('Tiempo (s)');
figure;plot(t,oxi_actual_red );
title('Lectura red (R)');
ylabel('Amplitud');
xlabel('Tiempo (s)');

%% Vectores x para graficar frecuencias
oxi_red=oxi_actual_red;
oxi_ir=oxi_actual_ir;
% Longitud de la TTF de los datos obtenidos
fourier_lenght=length(oxi_ir);
% Rango de frecuencia rads/s
w=0:2*pi/(fourier_lenght-1):2*pi;
% Frecuencia a Hertz
fre=(w*Fs)/(2*pi);


%% Dise�o del filtro passa bandas con ventana triangular para eliminar rizo
%48 ppm 
wp = 0.8/ (Fs/2);
%240 ppm
wr= 4/(Fs/2);
% Filtro pasa bandas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bandpass = fir1(100,[wp wr], 'bandpass', triang(101));
freqz(bandpass); title('filtro pasa bandas');
%% Eliminar componente directa de los datos num�ricos de R, IR
ir = oxi_ir - mean(oxi_ir);
red = oxi_red - mean(oxi_red);

%% Graficar la transformada de Fourier de las se�ales sin filtro ni DC
FFT_IR_ = fft(ir);
%graficas de las TFF
figure; plot(fre,abs(FFT_IR_));
ylabel('Magnitud');
xlabel('Frecuencia (Hz)');
title('Espectro Frecuencias infrarojo sin componente DC(IR)');

%%FFT for RED signal
FFT_RED_ = fft(red);
%graficas de las TFF
figure;plot(fre,abs(FFT_RED_));
ylabel('Magnitud');
xlabel('Frecuencia (Hz)');
title('Espectro Frecuencias rojo sin componente DC (R)');

%% Apllicar el filtro pasa bandas a las nuevas se�ales sin componente 
%% directa. El filtro pasa bandas elimina las frecuencias que no son de 
%% nuestro inter�s, aquellas que son menor a los 60ppm y mayores a 230ppm
p_banda_red=filter(bandpass,1,red);
p_banda_ir=filter(bandpass,1,ir);


%% Graficar la transformada de Fourier de las se�ales sin componente 
%% directa y filtradas 

FFT_IR = fft(p_banda_ir);
%graficas de las TFF
figure; plot(fre,abs(FFT_IR));
ylabel('Magnitud');
xlabel('Frecuencia (Hz)');
title('Espectro Frecuencias infrarojo filtrado(IR)');

%%FFT for RED signal
FFT_RED = fft(p_banda_red);
%graficas de las TFF
figure;plot(fre,abs(FFT_RED));
ylabel('Magnitud');
xlabel('Frecuencia (Hz)');
title('Espectro Frecuencias rojo filtrado(R)');


%% Sacamos promedio de max y mins de la se�al para poder realizar los 
%% los calculos de saturaci�n de ox�geno a trav�s de las diferencias entre
%% las amplitudes de la se�al

%% C�lculo de saturaci�n de ox�geno en la sangre

%% Encontrar los valores m�ximos de las se�ales, R e IR
vals_max_oxi1_red = abs(findpeaks(p_banda_red));
average_max_oxi1_red = mean(vals_max_oxi1_red);

% Findpeaks encuentra los m�ximos relativos, para su uso es necesario 
% hacer un promedio de los picos maximos
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
%% y considerando que 1 minuto (Hz)= 60 pulsaciones (1 cada segundo)
HEART_RATE = mean([max_IR, max_RED])*60;

%% Calcular la concentraci�n de ox�geno en la sangre aplicando un logaritmo
%% a la se�al, esto se puede encontrar en documentaci�n de diferentes f
%% fuentes donde se aplica una formula: R = log(ACred) / log(ACir) una
%% relaci�n entre las componentes directas y de ac de los valores del oxim.
% Calculo de Sp02
R_RED1 =log10(AC_oxi1_red);

R_IR1=log10(AC_oxi1_ir);

R1=R_RED1/R_IR1;
%% El siguiente paso es realizar el calculo con f�rmula 
%% C�lculado a trav�s de la formula de la comunidad de matlab
%% SpO2 = 100 - 25 * R 
SpO2 = 110 - 25*R1;

%% Graficar los valores de heart rate y de saturaci�n de ox�geno en la 
%% sangre a trav�s de una gr�fica 

% Graficas de pulso cardiaco
label = categorical({'Pulso cardiaco'});
subplot(1,2,1);
stem(label, HEART_RATE);
title('Pulso card�aco');
ylabel('Pulsos por minuto');
ylim([0 220]);

% Grafica de saturaci�n de oxigeno
label2 = categorical({'Saturaci�n de Ox�geno en la sangre'});
subplot(1,2,2);
stem(label2, SpO2);
title('Saturaci�n de Ox�geno en la sangre');
ylabel('Porcentaje de saturaci�n');
ylim([0 110]);






