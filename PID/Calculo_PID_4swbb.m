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

H=Rs2/(Rs1+Rs2);
%%%%%%%%%%Más variables%%%%%%%%

P=Vo^2/Ro;
D=Vo/(Vo+Vg);
tdpwm=D*Ts; %Esto cambia
td=tdpwm+tcntrl;
Il=P/(Vg*D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Propiedades del regulador%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fc=1/(Ts*30);
PM=45; %en grados
PM=PM*pi/180; %pasar a radianes

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
ylim([-180 370])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%-------Calcular regulador-----------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplicar por la ganancia
Tvuz=series(H,Tvuz);

%Lo de cambiar a p es solo si usas el espacio de estados, sino se puede
%calcular directamente el regulador y depues cambiarlo a z.

%nueva frecuencia en p (el prewrapping)

fc_s=(2/Ts)*tan(2*pi*fc*(Ts/2))/(2*pi); %en Hz
wc_s=fc_s*2*pi;
%pasar de discreta  a continua
opts = d2cOptions('Method','tustin','PrewarpFrequency',wc_s); 
Gvup_t=d2c(Tvuz,'tustin',opts);
figure('Name','Funcion de transferencia con tustin')
bodef(Gvup_t)


%Obtener los datos de bode en z
[mag_z,fase_z] = bode(Tvuz,fc*2*pi); %mag en ganancia natural y fase en grados
[mag_z1,fase_z1] = bode(Tvuz,0.1); %mag en ganancia natural y fase en grados
if(fase_z1>160) 
    fase_z=fase_z-360;
end
fase_z = (pi/180)*fase_z;%pasar a radianes



%CALCULAR EL PID Gz=Kp+Ki/(1-z^-1)+Kd(1-z^-1)
%En p GPID_p=Kp+(Ki*(1+p/wp))/(Ts*p)+(Kd*Ts*p)/(1+p/wp)
%Es mas facil usar la forma multiplicacion (ver en el libro cual es)

%Primero calcular el PD, que es el que regula la fase y mgnitud en la fc

wp=2/Ts;

% PI zero and high-frequency gain

wpd=wc_s/tan(PM-pi-fase_z+atan(wc_s/wp));

if(wp<wpd || wp<0 || wpd<0 )
    fprintf('Error al calcular PID\n');
    fprintf('El margen de fase debe estar entre %.4g y %.4g o varia la frecuancia de corte',(fase_z+pi)*180/pi , (fase_z+pi+pi/2-atan(wc_s/wp))*180/pi )
    return
end
fprintf('El margen de fase debe estar entre %.4g y %.4g',(fase_z+pi)*180/pi , (fase_z+pi+pi/2-atan(wc_s/wp))*180/pi )
Gpd = 1/mag_z*(sqrt(1+(wc_s/wp)^2))/(sqrt(1+(wc_s/wpd)^2));

%CALCULAR EL PI
wpi=1/20*wc_s;
Gpi=1;
% Calcular ganancias
Kp = Gpi*Gpd*(1+wpi/wpd-2*wpi/wp);
Ki = 2*Gpi*Gpd*wpi/wp;
Kd = Gpi*Gpd/2*(1-wpi/wp)*(wp/wpd-1);
% Funcion d transferencia
z = tf('z',Ts);
Gcz = Kp + Ki/(1-z^-1) + Kd*(1-z^-1);
figure('Name','TF de PID')
bodef(Gcz)
figure('Name','TF de todo con regulador')

T=series(Gcz,Tvuz);
bodef(T)
%No hace falta convertir el PID al espacio z, pq ya se ha hecho con una
%relacion, para mas informacion ver el libro 




function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end

