clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%
%%%%Parametros%%%%
%%%%%%%%%%%%%%%%%%

Ci=470*10^-6; 
Rci=62*10^-3; 
Rg=0;  
L=220*10^-6;
Rl=439*10^-3;  
Co=300*10^-6;
Rco=12.5*10^-3;
Vo=12;
Vg=4;
Ro= 100;

Ts=1/(40*10^3);

tcntrl=4e-6; %el retraso de ADC y regulador
Rs1=620*10^3;
Rs2=38.3*10^3;

H=Rs2/(Rs1+Rs2); % T da igual pq no lo vas a usar

%ADC
V_max_ADC=1;
V_min_ADC=0;
B_ADC=8;
q_adc=(V_max_ADC-V_min_ADC)/2^B_ADC;
%DPWM
V_max_DPWM=1;
B_DPWM=11;
q_pwm=((V_max_DPWM)/2^B_DPWM);

%margen de seguridad (debe de ser menor de 1)
alfa=0.7;

%margen de ganancia del sistema en db?
Gdb=8.52;


%%%%%%%%%%Más variables%%%%%%%%

P=Vo^2/Ro;
D=Vo/(Vo+Vg);
tdpwm=D*Ts; %Esto cambia
td=tdpwm+tcntrl;
Il=P/(Vg*D);



if( Gdb > 4.2-20*log10(alfa))
    fprintf("Se cumple B1, por lo cual puede que no haya ciclos limites");
else
    fprintf("No se cumple B1, por lo cual puede que haya ciclos limites");
end
