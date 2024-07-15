%Sirve para comprobar que la funcion de transferencia cuando en s y en z
%son las mismas, es decir la aproximacion se puede realizar.
clear ;
%close;
clc;
%%Parametros
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



%%%%%%%%%%%%%%Más variables%%%%%%%%
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
e1=[0;0];
e0=[0;0];
A=A1*D+A0*(1-D);
B=b1*D+b0*(1-D);
C=c1*D+c0*(1-D);
Ai = A^-1;
X=-Ai*B*[Vg];
Y=[-C*Ai*B+e1]*Vg;
F=(A1*X+b1*Vg)-A0*X+b0*Vg;
G=(c1*X+e1*Vg)-(c0*X+e0*Vg);
s = tf('s')
W=C*(s*eye(3)-A)^-1*F+G;
Gvd=[0 1]*W;
Gid=[1 0]*W;

%Promediado con resistencais parasitas (No se si esta bien)
s = tf('s');
A_p=(((1-D)*(Vg-Rg*Il+Vo))/(L*s+(Rg+Rl)*D+Rl*(1-D)))-Il;
B_p=(1-D)^2/(L*s+(Rg+Rl)*D+Rl*(1-D));
C_p=D*(1-D)/(L*s+(Rg+Rl)*D+Rl*(1-D));
Z_p=(Ro*Rco*Co*s+Ro)/(Co*s*(Ro+Rco)+1);
H1=feedback(Z_p,B_p);
Tu_par=series(H1,A_p);
figure('Name','Diferencia de SS y de promediado sin retardo')
hold on;
bode(Tu_par);
bode(Gvd)

%Comparar con la simulacion
%Abrir csv de simulink para ver las diferencias (sin retardos)
[mag1,phase1,wout1] = bode(Tu_par);
[mag2,phase2,wout2] = bode(Gvd);
fout=wout1/(2*pi);
fout2=wout2/(2*pi);
%leer el archivo
data=readmatrix('Bode_4sw_real.csv');
%crear el bode
f_psim=data(:,1);  
gain_psim=data(:,2);
phase_psim=data(:,3);
%Representacion y comparacion de bodes 
figure('Name','Diagrama de bode de Psim y Matlab')
subplot(2,1,1);
semilogx(fout,20*log10(mag1(:)),f_psim,gain_psim,fout2,20*log10(mag2(:)))
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
legend('Bode teórico','Bode Psim')
grid on
axis tight
subplot(2,1,2)
semilogx(fout,phase1(:),f_psim,phase_psim+360,fout2,phase2(:))
xlabel('Frecuencia (Hz)')
ylabel('Fase (grados)')
legend('Bode teórico con retardo','Bode Psim','Bode teórico sin retardo')
grid on
axis tight