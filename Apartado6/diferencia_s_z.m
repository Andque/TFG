clc;
clear
close all;
%Caracteristicas
C=200e-6;
rC=0.8e-0;
L=1e-6;
rL=30e-3;
Vg=5;
Vo=1.8;
Io=5;
Ts=1e-6;
D=Vo/Vg;
Dprime=1-D;
R=Vo/Io;
td=760e-9;
%METODO DE SS
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
figure("Name","BODE de SS")
bodef(Gvuz)
ylim([-370 10])

%%modelado promediado
Gvus=tf([Vg],[L*C L/R 1]);
figure("Name","BODE de promediado sin parasitas")
bodef(Gvus);
ylim([-370 10])

%Con las resisntencias parasitas
Gvusp=tf([rC*C*Vg Vg],[L*C (rC+rL)*C 1]);
figure("Name","BODE de promediado con parasitas")
bodef(Gvusp);
ylim([-370 10])

%Con la aproximacion de e^-tds
figure("Name","BODE de promediado sin parasitas y retardo ")
G_delay=tf(1,1,'InputDelay',td);
Tu=series(Gvus,G_delay);
h=bodeplot(Tu,'r--');
p = getoptions(h);
p.FreqUnits = 'Hz';
setoptions(h,p);
ylim([-370 10])
%Con la aproximacion de e^-tds
figure("Name","BODE de promediado con parasitas y retardo ")
G_delay=tf(1,1,'InputDelay',td);
Tup=series(Gvusp,G_delay);
h=bodeplot(Tup,'r--');
p = getoptions(h);
p.FreqUnits = 'Hz';
setoptions(h,p);
ylim([-370 10])

%Todo junto
figure("Name","Todo junto")
hold on

bode(Gvuz);
bode(Gvus);
bode(Gvusp);
bode(Tup);
bode(Tu);
ylim([-370 10])
%%Solo lo de digital
figure("Name","BODE de SS y promediado con parasitas y retardo")
hold on
bode(Tup);
bode(Gvuz)
ylim([-370 10])
%%Solo lo de digital
figure("Name","BODE de SS y promediado sin parasitas y retardo")
hold on
bode(Tu);
bode(Gvuz)
ylim([-370 10])


function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end