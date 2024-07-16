clc;
clear all;
close all;
Kp=3.09;
Ki=74.52e-3;
Kd=23.8;
Ts=1e-6;
z = tf('z',Ts);
Gcz= Kp + Ki/(1-z^-1) + Kd*(1-z^-1);
%para ver la estructura directa es abrir el Gcz y ver numerador y
%denominador y poner los numeradores a la izq y los den negados a la dch
%hay que ver que el primer valor del denominador sea 1, si no es 1 hay que
%multiplicar por (1/valor)/(1/valor).



%cascada, se puede usar sacando los polos y ceros de la Tf y despues lo
%ordenas como quieras.Importante poner la ganancia y todos los valores con
%signo contrario con el que aparecen y ponerlo como -1/valor.
[numGc,denGc] = tfdata(Gcz);
numGc = cell2mat(numGc); %combiertes todo a una misma matriz
denGc = cell2mat(denGc);
[zeros,polos,ganancia] = tf2zp(numGc,denGc);
