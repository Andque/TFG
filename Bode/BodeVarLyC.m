D=0.5;
Vin=10;
Vout=10;
R=20;
L=1;
Il=1;
C=1 ;
A=(1-D)*(Vin+Vout)*R;
B=R*L*Il;
F=R*((1-D)^2);
D=C*R*L;
E=L;
H= tf([-B A],[D E F]);
%%%%%%%
D1=0.5;
Vin1=10;
Vout1=10;
R1=20;
L1=5.625*(10^(-3));
Il1=0.4;
C1=40*(10^(-6)) ;
A1=(1-D1)*(Vin1+Vout1)*R1;
B1=R1*L1*Il1;
F1=R1*((1-D1)^2);
D1=C1*R1*L1;
E1=L1;
H1= tf([-B1 A1],[D1 E1 F1]);
%%%%%%%
D2=0.5;
Vin2=10;
Vout2=10;
R2=20;
L2=5.625*(10^(-12));
Il2=0.4;
C2=40*(10^(-16)) ;
A2=(1-D1)*(Vin1+Vout1)*R1;
B2=R1*L2*Il2;
F2=R2*((1-D2)^2);
D2=C2*R1*L2;
E2=L2;
H2= tf([-B2 A2],[D2 E2 F2]);
bode(H,H1,H2);