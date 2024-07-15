close all;
clear;
clc;
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
fsw=1/Ts;

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

fc=fsw/60; %aproximado
PM_g=70; %en grados
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
hold on
bodef(Tvuz)
Tvuz=series(Tvuz,H);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------Calcular regulador-----------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Paso 1:Pasamos la funcion de transferencia del convertidor de z a p

fc_p=(2/Ts)*tan(2*pi*fc*(Ts/2))/(2*pi); %en Hz %Pasamos la frecuencia angular de corte a la frecuencia angular en el dominio p 
wc_p=fc_p*2*pi;
opts = d2cOptions('Method','tustin','PrewarpFrequency',wc_p); 
Gvup_t=d2c(Tvuz,'tustin',opts);
figure("Name","Funcion de transferencia pasada a p con tustin")
bode(Gvup_t);

%Esto nos da la ganancia en unidades absolutas y la fase en grados, para la frecuencia de corte y para 0 Hz
[GTup,FTup]=bode(Gvup_t,fc_p*2*pi); 
[GTup0,FTup0]=bode(Gvup_t,0.1); 

%Para cambiar la fase a una negativa (es con la que se trabaja)
if(FTup0>200) %es decir empieza por un numero mayor que 180, lo normal es que empiece en 0
    FTup=FTup-360;
end

%Usar el metodod de la K para calcular el regulador
K=tan(1/4*(pi/2+PM_r-FTup*(pi/180)));

if(K<1)
disp("Error para calcular el tipo 3");
return;
end

R11=20*10^3; %Valor predeterminado, se puede cambiar
R1=R11/(K^2-1);
C1_cond=1/(wc_p*R1*K);
C23_cond=GTup/(wc_p*R11)*K^2;
C3_cond=C23_cond/(K^2);
C2_cond=C23_cond-C3_cond;
R2=K/(wc_p*C2_cond);
fz=1/(2*pi*R2*C2_cond);
fp=1/(2*pi*R2*((C2_cond*C3_cond)/C23_cond));
display("En el tipo 3 con K " +K + "y fc de " + wc_p/(2*pi)+": la frecuencia de polo es " +fp+ " y la de zero " +fz);

%Ver la funcion de transferencia en p del tipo 3
p = tf('p');
Gcp3=( (p*C1_cond*(R11+R1)+1) * (p*C2_cond*R2+1) ) / ( (p*R11*(C2_cond+C3_cond)) * (p*C1_cond*R1+1) * (p*R2*((C2_cond*C3_cond)/(C2_cond+C3_cond))+1) );
figure('Name',"FT de sistema en p")
hold on;
Gvup_tot=series(Gcp3,Gvup_t);
bode(Gcp3,Gvup_t,Gvup_tot);

%Transformar del dominio p a dominio z con el mapeado bilinear
opt=c2dOptions('Method','tustin','PrewarpFrequency',wc_p);
Gcz3 = c2d(Gcp3,Ts,opt);

%Ver los ceros y polos de este sistema (Para poder hacer la forma paralela del regulador)
[numGc,denGc] = tfdata(Gcz3); %sacas el denominador y numerador de la ecuacion de trasnferencia
numGc = cell2mat(numGc); %combiertes todo a una misma matriz
denGc = cell2mat(denGc);
disp("La FT del tipo 3 tiene la forma de: ") %Forma directa
disp("numerador:");
disp(numGc);
disp("denominador:");
disp(denGc)

[numGc,denGc] = eqtflength(numGc,denGc); %Pones todo con la misma longuitud

% Ver la forma paralela 
% en plan, seria separar todo en sus fracciones parciales
% y ya, excepto si alguna tiene numeros complejos, que lo tendrias que
% multiplicar entre si para quitarlo.
[rparalel,poleparalel,kparalel] = residuez(numGc,denGc); %Lo usas para calcular la fraccion en sumas
disp("ecuacion en fraciones parciales:");
disp("coef de arriba:");
disp(rparalel);
disp("polos :");
disp(poleparalel);
disp("ganancia:");
disp(kparalel);

%%%
figure('NAME',"TF de regulador, planta sin regular y regulada");
hold on;
bode(Tvuz);
bode(Gcz3);
Tz3=series(Tvuz,Gcz3);
bode(Tz3);
[Gm,Pm,Wcg,Wcp] = margin(Tz3); %Usado para saber si es 

if(Pm<=0)
    disp("El sistema no es estable, elige otro margen de fase o Fc");
    return;
end

%Estructura paralela para el FPGA 
%Vamos a poner la misma notacion que en el dibujo del word
%Esto solo vale para nuestro controlador, ya que puede que tenga otras formas
K1=real(rparalel(1));
K2=real(rparalel(2));
K3=real(rparalel(3));
C1=real(poleparalel(1));
C2=real(poleparalel(2));
C3=real(poleparalel(3));
K4=kparalel;
figure('NAME',"regulador // y directo");
[Gcz3_num,Gcz3_den]= residuez(rparalel,poleparalel,kparalel);
hold on;
bode( Gcz3)

z = tf('z',Ts);
Gcz3_p=(K1/(1-C1*z^-1))+(K2/(1-C2*z^-1))+K4 + (K3/((1-C3*z^-1)^2));
bode(Gcz3_p);

% %%%ABRIR CSV DE BODE, se hizo con H=1;
% 
% [mag1,phase1,wout1] = bode(Gcz3_p);
% fout=wout1/(2*pi);

% %leer el archivo, esto se hizo sin contar la H, asi que para que salga
% %igual H=1
% data=readmatrix('Bode_tipo3_prueba.csv');
% %crear el bode
% f_psim=data(:,1);  
% gain_psim=data(:,2);
% phase_psim=data(:,3);
% %Representacion y comparacion de bodes 
% figure('Name','Diagrama de bode de Psim y Matlab')
% subplot(2,1,1);
% semilogx(fout,20*log10(mag1(:)),f_psim,gain_psim)
% xlabel('Frecuencia (Hz)')
% ylabel('Magnitud (dB)')
% legend('Bode teórico','Bode Psim')
% grid on
% axis tight
% subplot(2,1,2)
% semilogx(fout,phase1(:),f_psim,phase_psim)
% xlabel('Frecuencia (Hz)')
% ylabel('Fase (grados)')
% legend('Bode teórico','Bode Psim')
% grid on
% axis tight
% 
% 


function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end

