clc;
clear;
close;
%sirve para el que estamos analizando.Para un bick seia distinto la
%ecuacion del qvo_DPWM
%variables
Vfs=2; %valor maximo al que llega el ADC, el minimo se supone que es 0
B_ADC=8;
B_PWM=10;
Vo=1.8; %Voltaje de salida
Vg=5;   %voltaje de entrada
D=Vo/Vg; %D
Nr=1;   %Valor maximo al que llega el DPWM %Por facilidad vamos a suponer que es siempre 1, despues en el regulador se solucionara esto
H=1;    %TF del sensado


%calculos
qu=Nr/2^B_PWM; %La resolucion de el DPWM
qvo_ADC=Vfs/(2^B_ADC*H);
D_cuant=Qn(D,qu,Nr);
V_ref_cuant=Qn(Vo*H,qvo_ADC,Vfs);

qvo_DPWM=qu/Nr*Vg;
display("Qdpwm: "+qvo_DPWM+"   Qadc: "+qvo_ADC);
if (qvo_DPWM>=qvo_ADC)
    fprintf("Hay que subir la resolucion al DPWM");
end
if(qvo_DPWM<qvo_ADC) 
    fprintf("Es valida la resolucion al DPWM");
end

%%Generar una funcion que te ponga graficamente los VD y los VADC
q_adc=linspace(0,Vfs-(Vfs/2^B_ADC),2^B_ADC);
figure('Name','Valores de qADC y qDPWM')
hold on;
xline(q_adc,'m','LineWidth',1)
xline(V_ref_cuant,'r','LineWidth',2)

%la funcion de qvo_dpwm
q_dpwm=linspace(0,Vg,2^B_ADC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
Vs_dpwm = zeros(1,2^B_PWM);
n=1;
D_cuant_n=0;
while n<=2^B_PWM
 Vs_dpwm(1,n)=D_cuant_n*Vg;
 n=n+1;
  D_cuant_n=D_cuant_n+qu;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



xline(Vs_dpwm,'k','LineWidth',1)
xlabel('Vs, pink=qvs\_ADC, black=qvs\_DPWM');
%%calcular el V de referencia cuantificado
function n = Qn(Dx,qd,Vfs)
n=0;
if(n+qd>Dx)
    return
end
n=qd;
while(n+qd<=Dx && n<(Vfs-qd))
    n=n+qd;
end
return;
end