clc;
clear
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
td=760e-9;
%Propiedades del regulador
fc=1/(Ts*10);
PM=45; %en grados
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

%nueva frecuencia de corte
fc_s=(2/Ts)*tan(2*pi*fc*(Ts/2))/(2*pi); %en Hz
fp=2/(Ts);
%pasar de discreta  a continua
opts = d2cOptions('Method','tustin','PrewarpFrequency',fp); 
Gvup_t=d2c(Gvuz,'Tustin',opts); 
Gvup_d=d2c(Gvuz,'ZOH'); %Mal pq no mantiene la estabilidad
%Obtener los datos de la TF
%obtener el numerador y denominador separados
[numz,denz]=tfdata(Gvuz);
numz2=zeros(1,size(numz{1,1},2));
for i=1:size(numz{1,1},2)
    numz2(1,i)=numz{1,1}(1,i);
end

denz2=zeros(1,size(denz{1,1},2));
for i=1:size(denz{1,1},2)
    denz2(1,i)=denz{1,1}(1,i);
end
%crear una funcion fantasma con la tf !!!!!!!!!!!!!
%!!!!!!!!!!!!!!!HACERLO A MANO, osea con el mmismo grado que tienes
f=@(x) (numz2(1,1)*x^2+numz2(1,2)*x^1+numz2(1,3)*x^0)/...
    (denz2(1,1)*x^2+denz2(1,2)*x^1+denz2(1,3)*x^0);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
s=tf('s');
%Probar a hacer el ZOH a mano
figure('Name','Matlab vs a mano de Euler Delante')
hold on
bode(Gvup_d)
Gvup_dm=f(1+s*Ts);
h = bodeplot(Gvup_dm);
setoptions(h,'PhaseMatching','on','PhaseMatchingFreq',1,'PhaseMatchingValue',0,'FreqUnits','Hz','Xlim',[1000 ,1/Ts]);


figure('Name','Matlab vs a mano de Tustin')
hold on
bode(Gvup_t)
Gvup_tm=f((2+s*Ts)/(2-s*Ts));
h = bodeplot(Gvup_tm);
setoptions(h,'PhaseMatching','on','PhaseMatchingFreq',1,'PhaseMatchingValue',0,'FreqUnits','Hz','Xlim',[1000 ,1/Ts]);

figure('Name','Matlab  a mano de Euler detras')
Gvup_am=f(1/(1-s*Ts));
bode(Gvup_am)

figure('Name','Euler y Tustin comparacion de funcion de matlab')
hold on
bode(Gvup_d);
h = bodeplot(Gvup_t);
setoptions(h,'PhaseMatching','on','PhaseMatchingFreq',1,'PhaseMatchingValue',0,'FreqUnits','Hz','Xlim',[1000 ,1/Ts]);

figure('Name','Euler y Tustin comparacion de funcion a mano')
hold on
bode(Gvup_am);
bode(Gvup_tm);
h = bodeplot(Gvup_dm);
setoptions(h,'PhaseMatching','on','PhaseMatchingFreq',1,'PhaseMatchingValue',0,'FreqUnits','Hz','Xlim',[1000 ,1/Ts]);

[mag,phase] = bode(Gvup_t,fc_s*2*pi);
phase=phase-360;
mag=20*log10(mag);



function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end

