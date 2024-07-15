%Sirve para comprobar que la funcion de transferencia cuando en s y en z
%son las mismas, es decir la aproximacion se puede realizar.
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

%Sin retardo espacio de Estado

% Poner las matrices A1 es cuando el Il sube y A0 cuando baja.
A1_s=[ -(1/(Ci*(Rg+Rci))) 0 -(Rg/(Ci*(Rg+Rci))); 0 -(1/(Co*(Ro+Rco))) 0; Rg/(L*(Rg+Rci)) 0 -((Rg*Rl+Rci*Rg+Rl*Rci)/(L*(Rg+Rci))) ];
A0_s=[ -(1/(Ci*(Rg+Rci))) 0 0; 0 -(1/(Co*(Ro+Rco))) Ro/(Co*(Ro+Rco)); 0 -(Ro/(L*(Ro+Rco))) -((Rco*Rl+Rco*Ro+Rl*Ro)/(L*(Ro+Rco)))];
b1_s=[ 1/(Ci*(Rg+Rci)); 0; Rci/(L*(Rg+Rci))];
b0_s=[ 1/(Ci*(Rg+Rci)); 0; 0];
c1_s=[ 0 0 1; 0 Ro/(Rco+Ro) 0];
c0_s=[ 0 0 1; 0 (Ro/(Ro+Rco)) (Ro*Rco/(Ro+Rco))];
e1_s=[0;0];
e0_s=[0;0];
A_s=A1_s*D+A0_s*(1-D);
B_s=b1_s*D+b0_s*(1-D);
C_s=c1_s*D+c0_s*(1-D);
Ai_s = A_s^-1;
X_s=-Ai_s*B_s*(Vg);
Y_s=(-C_s*Ai_s*B_s+e1_s)*Vg;
F_s=(A1_s*X_s+b1_s*Vg)-A0_s*X_s+b0_s*Vg;
G_s=(c1_s*X_s+e1_s*Vg)-(c0_s*X_s+e0_s*Vg);
s = tf('s');
W_s=C_s*(s*eye(3)-A_s)^-1*F_s+G_s;
Gvd_s=[0 1]*W_s;
Gid_s=[1 0]*W_s;

%%En analogico con el retardo sin resistencias parasitas!!! habria que
%%calcular con resistencias parasitas para ver si es lo mismo, y no lo es
%%debido a que las resistencias parasitas te bajan la fase.
A=(1-D)*(Vg+Vo)*Ro;
B=Ro*L*Il;
F=Ro*((1-D)^2);
De=Co*Ro*L;
E=L;  
Tu= tf([-B A],[De E F]); %solo del convertidor
G_delay=tf(1,1,'InputDelay',td);
Tu_re=series(Tu,G_delay);


%En analogico con resitencias parasitas
s = tf('s');
A_p=(((1-D)*(Vg-Rg*Il+Vo))/(L*s+(Rg+Rl)*D+Rl*(1-D)))-Il;
B_p=(1-D)^2/(L*s+(Rg+Rl)*D+Rl*(1-D));
C_p=D*(1-D)/(L*s+(Rg+Rl)*D+Rl*(1-D));
Z_p=(Ro*Rco*Co*s+Ro)/(Co*s*(Ro+Rco)+1);
H=feedback(Z_p,B_p);
Tu_par=series(H,A_p);





%%Plot el csv

Bode_psim=readmatrix('Bode_sin_res.csv');

%crear una matriz para comparar el bode con csv sin resistencias
%parasitarias
[m,n]=size(Bode_psim);
X=zeros(m,n);
for i=1:m
[mag,phase,wout] = bode(Tu,Bode_psim(i,1)); 
X(i,1)=wout/(2*pi);
X(i,2)=20*log10(mag);
X(i,3)=phase-360;
end

%%
figure('Name','Diferencia de promediado vs psim sin parasitas')
hold on
tiledlayout(2,1)
% Top plot
nexttile
semilogx(Bode_psim(:,1),Bode_psim(:,2))
hold on
semilogx(X(:,1),X(:,2))
hold on
title('Magnitud')

% Bottom plot
nexttile
semilogx(Bode_psim(:,1),Bode_psim(:,3))
hold on
semilogx(X(:,1),X(:,3))
hold on
title('Fase')



%crear una matriz para comparar el bode con csv con resistencias
Bode_psim_res=readmatrix('Bode_con_res.csv');

%parasitarias
[m_res,n_res]=size(Bode_psim_res);
X_res=zeros(m_res,n_res);
for i=1:m_res
[mag_res,phase_res,wout_res] = bode(Tu_par,Bode_psim_res(i,1));  %Tu_par para promediado y Gv_s para SS
X_res(i,1)=wout_res/(2*pi);
X_res(i,2)=20*log10(mag_res);
X_res(i,3)=phase_res-360;
end
%%
figure('Name','Diferencia de promediado vs psim sin retardo y con parasitas')
hold on
tiledlayout(2,1)
% Top plot
nexttile
semilogx(Bode_psim_res(:,1),Bode_psim_res(:,2))
hold on
semilogx(X_res(:,1),X_res(:,2))
hold on
title('Magnitud')

% Bottom plot
nexttile
semilogx(Bode_psim_res(:,1),Bode_psim_res(:,3))
hold on
semilogx(X_res(:,1),X_res(:,3))
hold on
title('Fase')

