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
%%%%%%%%%%Más variables%%%%%%%%

P=Vo^2/Ro;
D=Vo/(Vo+Vg);
tdpwm=D*Ts; %Esto cambia
td=tdpwm+tcntrl;
Il=P/(Vg*D);
H=Rs2/(Rs1+Rs2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Propiedades del regulador%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc=372; %aproximado
PM_g=40; %en grados
PM_r=PM_g*pi/180; %pasar a radianes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Funcion de transf en la peor condicion%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Espacio de Estado

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
Xdown = ((eye(3)-expm(A1*D*Ts)*expm(A0*Dprime*Ts))^-1)*(-expm(A1*D*Ts)*A0i*(eye(3)-expm(A0*Dprime*Ts))*b0+ -A1i*(eye(3)-expm(A1*D*Ts))*b1)*[Vg];
Phi = expm(A0*(Ts-td))*expm(A1*D*Ts)*expm(A0*(td-D*Ts));
%Asumes que el Nr es 1 (voltage del DPWM) y H=1, esto no es del todo cierto
%pero en el escalado 
%del regulador lo vas a  solucionar se va a dividir entre 1/Nr para que si sea cierto (Para mayor comodidad)
gamma = expm(A0*(Ts-td))*((A1-A0)*Xdown + (b1-b0)*[Vg])*Ts;
delta = c1; %Aqui pone la c0 si el muestreado ocurre durante el periodo 0, o c1 si el muestreado ocurre durante el periodo 1
%Despues extraer la funcion que queremos
sys = ss(Phi,gamma,delta(2,:),0,Ts); % Construir un objeto que sea del "State Space"
Tvuz = tf(sys); %Obtener la funcion de transferencia de dicho objeto
figure('Name','Funcion de transferencia del convertidor H=1 y NR=1')
bodef(Tvuz)
figure('Name','Funcion de transferencia del convertidor NR=1')
Tvuz=series(Tvuz,H);
bodef(Tvuz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------Calcular regulador-----------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lo de cambiar a p es solo si usas el espacio de estados, sino se puede
%calcular directamente el regulador y depues cambiarlo a z.

%nueva frecuencia en p (el prewrapping)

fc_s=(2/Ts)*tan(2*pi*fc*(Ts/2))/(2*pi); %en Hz
wc_s=fc_s*2*pi;
wp=2/Ts;
%pasar de discreta  a continua
opts = d2cOptions('Method','tustin','PrewarpFrequency',wc_s); 
Gvup_t=d2c(Tvuz,opts);
figure("Name","Funcion de transferencia pasada a p con tustin y en z")
bode(Gvup_t,Tvuz);
%Obtener los datos de bode en z
[mag,phase,wout] = bode(Tvuz);
[mag_z1,fase_z1] = bode(Tvuz,0.1); %mag en ganancia natural y fase en grados
mag = squeeze(mag);
phase = squeeze(phase);
if(fase_z1>160) 
    phase=phase-360;
end
freq_deseada = interp1(phase, wout, (180+PM_g)-360)/(2*pi);               % Find Desired Frequency
[mag_z2,fase_z2] = bode(Gvup_t,freq_deseada*2*pi); %mag en ganancia natural y fase en grados

Kc=1/mag_z2;

%calcular el valor de la integral que no debe de interferir con la el
%la frecuencia de corte, por lo cual vamos a decir que 
Ti=0.00001*freq_deseada*2*pi;
s = tf('s')
Gi=(s*Ti+1)/s;
Gki=series(series(Gi,Kc),1/Ti);
Kc=Kc/Ti;
figure("Name","FT de regulador")
bode(Gki)

Gvup_t_ki=series(Gvup_t,Gki);
figure("Name","Todo junto")
bode(Gvup_t_ki)

%Calcular los valores de Kp y Ki en p
Kp=Kc*Ti; 
Ki=Kc;

%Caclular los valores de Kp y Ki en z:

opts = c2dOptions('Method','tustin','PrewarpFrequency',wc_s);
Gc_z = c2d(Gki,Ts,opts);

%Ver los ceros y polos de este sistema (Para poder hacer la forma paralela del regulador)
[numGc,denGc] = tfdata(Gc_z); %sacas el denominador y numerador de la ecuacion de trasnferencia
numGc = cell2mat(numGc); %combiertes todo a una misma matriz
denGc = cell2mat(denGc);
disp("La FT del tipo PI tiene la forma de: ") %Forma directa
disp("numerador:");
disp(numGc);
disp("denominador:");
disp(denGc)


% Ver la forma paralela 
% en plan, seria separar todo en sus fracciones parciales
[PI_num,PI_den,PI_K] = residuez(numGc,denGc); %Lo usas para calcular la fraccion en sumas

Kp_z1=-numGc(1,2);
Ki_z1=numGc(1,1)-Kp_z1;

figure("Name","Funcion de transferencia PI en z")

bode(Gc_z)

function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end
