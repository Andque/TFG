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

%frecuencia donde el valor de la fase de la planta y regulador es -180
fx=1.48e3;
%pasarla a rad/s
fx=fx*2*pi;
%%%%%%%%%%Más variables%%%%%%%%

P=Vo^2/Ro;
D=Vo/(Vo+Vg);
tdpwm=D*Ts; %Esto cambia
td=tdpwm+tcntrl;
Il=P/(Vg*D);

%%Espacio de Estado con retardo

% Poner las matrices A1 es cuando el Il sube y A0 cuando baja.
A1=[ -(1/(Ci*(Rg+Rci))) 0 -(Rg/(Ci*(Rg+Rci))); 0 -(1/(Co*(Ro+Rco))) 0; Rg/(L*(Rg+Rci)) 0 -((Rg*Rl+Rci*Rg+Rl*Rci)/(L*(Rg+Rci))) ];
A0=[ -(1/(Ci*(Rg+Rci))) 0 0; 0 -(1/(Co*(Ro+Rco))) Ro/(Co*(Ro+Rco)); 0 -(Ro/(L*(Ro+Rco))) -((Rco*Rl+Rco*Ro+Rl*Ro)/(L*(Ro+Rco)))];
b1=[ 1/(Ci*(Rg+Rci)); 0; Rci/(L*(Rg+Rci))];
b0=[ 1/(Ci*(Rg+Rci)); 0; 0];
c1=[ 0 0 1; 0 Ro/(Rco+Ro) 0];
c0=[ 0 0 1; 0 (Ro/(Ro+Rco)) (Ro*Rco/(Ro+Rco))];
%Calculos si el DPWM es de trailing edge, la matriz A1 y A0 debe de ser invertible
Dprime=1-D;
A1i = A1^-1;
A0i = A0^-1;
Xdown = ((eye(3)-expm(A1*D*Ts)*expm(A0*Dprime*Ts))^-1)*(-expm(A1*D*Ts)*A0i*(eye(3)-expm(A0*Dprime*Ts))*b0+ -A1i*(eye(3)-expm(A1*D*Ts))*b1)*(Vg);
Phi = expm(A0*(Ts-td))*expm(A1*D*Ts)*expm(A0*(td-D*Ts));
%Asumes que el Nr es 1 (voltage del DPWM) y H=1, esto no es del todo cierto
%pero en el escalado 
%del regulador lo vas a  solucionar se va a dividir entre 1/Nr para que si sea cierto (Para mayor comodidad)
gamma = expm(A0*(Ts-td))*((A1-A0)*Xdown + (b1-b0)*(Vg))*Ts;
delta = c1; %Aqui pone la c0 si el muestreado ocurre durante el periodo 0, o c1 si el muestreado ocurre durante el periodo 1
%Despues extraer la funcion que queremos
sys = ss(Phi,gamma,delta(2,:),0,Ts); % Construir un objeto que sea del "State Space"
Tvuz = tf(sys); %Obtener la funcion de transferencia de dicho objeto
Tvuz=series(Tvuz,H);
bode(Tvuz)

%Ver cual es la ganancia cuando el tiene 180º

[mag1,phase1,wout1] = bode(Tvuz,fx);

if((4/pi) *mag1*q_pwm<q_adc*alfa)
    fprintf("Se cumple B1, por lo cual puede que no haya ciclos limites");
else
    fprintf("No se cumple B1, por lo cual puede que haya ciclos limites");
end
