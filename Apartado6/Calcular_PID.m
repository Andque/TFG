clc;
clear all;
close all;
%Cte de sistema
C=200e-6;
rC=0.8e-3;
L=1e-6;
rL=30e-3;
Vg=5;
Vo=1.8;
Io=5;
Ts=1e-6;
D=Vo/Vg;
Dprime=1-D;
R=Vo/Io;
tcntr=400e-9;
td=tcntr+D*Ts;
%Propiedades del regulador
fc=1/(Ts*10);
PM=45; %en grados
PM=PM*pi/180; %pasar a radianes
%calcular TF usando SS
A1 = [ -(rC+rL)/L -1/L; 1/C 0 ];
A0 = A1;
b1 = [ 1/L rC/L; 0 -1/C ];
b0 = [ 0 rC/L; 0 -1/C];
c1 = [1 0; rC 1];
c0 = c1;
A1i = A1^-1;
A0i = A0^-1;
Xdown = ((eye(2)-expm(A1*D*Ts)*expm(A0*Dprime*Ts))^-1)*...
(-expm(A1*D*Ts)*A0i*(eye(2)-expm(A0*Dprime*Ts))*b0+...
-A1i*(eye(2)-expm(A1*D*Ts))*b1)*[Vg;Io];
Phi = expm(A0*(Ts-td))*expm(A1*D*Ts)*expm(A0*(td-D*Ts));
gamma = expm(A0*(Ts-td))*((A1-A0)*Xdown + (b1-b0)*[Vg;Io])*Ts;
delta = c0;
sys = ss(Phi,gamma,delta(1,:),0,Ts);
Giuz = tf(sys);
sys = ss(Phi,gamma,delta(2,:),0,Ts);
Gvuz = tf(sys);
bodef(Gvuz)

%-------Calcular regulador-----------

%nueva frecuencia en p (el prewrapping)
fc_s=(2/Ts)*tan(2*pi*fc*(Ts/2))/(2*pi); %en Hz
wc_s=fc_s*2*pi;
%pasar de discreta  a continua
opts = d2cOptions('Method','tustin','PrewarpFrequency',wc_s); 
Gvup_t=d2c(Gvuz,'tustin',opts);

%Obtener los datos de bode en z
[mag_z,fase_z] = bode(Gvuz,fc*2*pi);
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
    fprintf('El margen de fase debe estar entre %.4g y %.4g',(fase_z+pi)*180/pi , (fase_z+pi+pi/2-atan(wc_s/wp))*180/pi )
    return
end
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

T=series(Gcz,Gvuz);
bodef(T)


function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end

