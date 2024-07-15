clc;
clear;
close;
%sirve para el que estamos analizando.Para un bick seia distinto la
%ecuacion del qvo_DPWM
%variables
Vfs=1; %valor maximo al que llega el ADC, el minimo se supone que es 0
B_ADC=8;
B_PWM=11;
Vo=12; %Voltaje de salida
Vg=12;   %voltaje de entrada maximo y minimo!!!!!!
D=Vo/(Vg+Vo); 
Nr=1;   %Valor maximo al que llega el DPWM %Por facilidad vamos a suponer que es siempre 1, despues en el regulador se solucionara esto
H=38.3/(38.3+620);    %TF del sensado


%calculos
qu=Nr/2^B_PWM; %La resolucion de el DPWM
qvo_ADC=Vfs/(H*(2^B_ADC));
D_cuant=Qn(D,qu,Nr);
V_ref_cuant=Qn(Vo,qvo_ADC,Vfs);

qvo_DPWM=Vg*qu/((1-D_cuant)^2 *Nr); 
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
D_dpwm = zeros(1,2^B_PWM);
Vs_dpwm = zeros(1,2^B_PWM);
n=2;
while n<=2^B_PWM
 D_dpwm(1,n)=qu*(n-1);
 n=n+1;
end

n=1;
D_cuant_n=0;
while n<=2^B_PWM
 Vs_dpwm(1,n)=(H*Vg*D_cuant_n)/(1-D_cuant_n);
 n=n+1;
  D_cuant_n=D_cuant_n+qu;
end

xline(Vs_dpwm,'k','LineWidth',1)
xlabel('Vs, pink=qvs\_DPWM, black=qvs\_ADC');
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