clc;
clear all;
close all;
K1=0.042909822280111;
K2=12.656362606545766;
K3=-3.160009652371582;
K4=-6.465264390256411;
C1=1;
C2=0.670560975909576;
C3=0.670560962239658;

syms x
Ki=limit( (x-1) *( (K1/(1-C1*x^-1))+(K2/(1-C2*x^-1))+K4 + (K3/((1-C3*x^-1)^2)) ), 1);
Ki_nice= double(Ki);