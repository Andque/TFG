clc
close all
clear all
%hay que ver en fc y en casi 0Hz la diferencia en % y dependiendo de el
%maximo permitido se escoge un numero de bits u otro.

%Valores constantes
V_max_ADC=1;
V_min_ADC=0;
B_ADC=8;
Nr=2047; %32767Poner el NR de verdad que en este paso se usa 
Kp=8.994696290258627;
Kd=1.25;
Ki=0.108832338837131;
Ts=2.5e-05;
fc= 1/(Ts*20); %frecuencia de corte

%%Especifico de este paso

B_min_r=3;  %minimo numero de bits para los valores del regulador MINIMO 2!!!, 3 si quieres asegurarte que la funcion funcione bien
B_max_r=20; %maximo numero de bits permitidos
Er_m_f0=10; %Porcentaje de diferencia en la magnitud en frecuancia 0
Er_f_fc=1;  %Porcentaje de diferencia en la fase en frecuancia fc
Er_m_fc=0.56;  %Porcentaje de diferencia en la magnitud en frecuancia fc

%PID en // con el escalado!! Que en este caso es multiplicar por Nr y qvs
qvs=(V_max_ADC-V_min_ADC)/(2^B_ADC);
Kp_e=Kp*Nr*qvs;
Kd_e=Kd*Nr*qvs;
Ki_e=Ki*Nr*qvs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Calculo de bits minimos   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Introducir la funcion de trasnferencia del PID sin cuantificar y con el
%escalado!!!!
z = tf('z',Ts);
Gcz_ncuant = Kp_e + Ki_e/(1-z^-1) + Kd_e*(1-z^-1);
figure('Name','TF de PID escalada y sin cuantificar')
bodef(Gcz_ncuant);
[mag_ncuant_0,phase_ncuant_0] = bode(Gcz_ncuant,0.1);
[mag_ncuant_fc,phase_ncuant_fc] = bode(Gcz_ncuant,fc*2*pi);

%calcular con los distintos valores de Kp,Ki y Kd cuantificados cual es el
%k hace que el error sea el maximo permitido con el minimo numero de bits.

%BUCLES!!!
for B_Kp = B_min_r:B_max_r
Kp_cuant=Qn(Kp_e,B_Kp); %Cuantificar Kp con los bits correspondientes
    for B_Kd = B_min_r:B_max_r
    Kd_cuant=Qn(Kd_e,B_Kd);%Cuantificar Kd con los bits correspondientes
        for B_Ki = B_min_r:B_max_r
         Ki_cuant=Qn(Ki_e,B_Ki);%Cuantificar Kd con los bits correspondientes
% Crear una nueva TF con los nuevos valores
         z = tf('z',Ts);
         Gcz_cuant = Kp_cuant.xq + Ki_cuant.xq /(1-z^-1) + Kd_cuant.xq *(1-z^-1);

         [mag_cuant_0,phase_cuant_0] = bode(Gcz_cuant,0.1);
         [mag_cuant_fc,phase_cuant_fc] = bode(Gcz_cuant,fc*2*pi);

%         Ver si se cumplen las 3 condiciones
         Er_m_f0_real=abs (abs(mag_ncuant_0-mag_cuant_0) / mag_ncuant_0) *100;
         Er_m_fc_real=abs( abs(mag_ncuant_fc-mag_cuant_fc) / mag_ncuant_fc * 100);
         Er_f_fc_real=abs( abs(phase_ncuant_fc-phase_cuant_fc) / phase_ncuant_fc * 100);
    %fprintf("\n Kp=  %g B_Kp= %u Forma: %s",Kp_cuant.xq,B_Kp,Kp_cuant.s);
    %fprintf("\n Kd=  %g B_Kd= %u Forma: %s",Kd_cuant.xq,B_Kd,Kd_cuant.s);
    %fprintf("\n Ki=  %g B_Ki= %u Forma: %s",Ki_cuant.xq,B_Ki,Ki_cuant.s);
    %fprintf("\n Los errores son:");
    %fprintf("\n Error magnitud en 0: %d",Er_m_f0_real);
    %fprintf("\n Error magnitud en fc: %d",Er_m_fc_real);
    %fprintf("\n Error fase en fc: %d",Er_f_fc_real);
         
         if( Er_m_f0_real<Er_m_f0 && Er_m_fc_real<Er_m_fc && Er_f_fc_real<Er_f_fc  )
            break
         end 
        end
    if( Er_m_f0_real<Er_m_f0 && Er_m_fc_real<Er_m_fc && Er_f_fc_real<Er_f_fc )
      break
    end
    end


if( Er_m_f0_real<Er_m_f0 && Er_m_fc_real<Er_m_fc && Er_f_fc_real<Er_f_fc )
  fprintf("\n El valor para que se cumpla los errores maximos permitidos son");
    fprintf("\n Kp=  %g B_Kp= %u Forma: %s",Kp_cuant.xq,B_Kp,Kp_cuant.s);
    fprintf("\n Kd=  %g B_Kd= %u Forma: %s",Kd_cuant.xq,B_Kd,Kd_cuant.s);
    fprintf("\n Ki=  %g B_Ki= %u Forma: %s",Ki_cuant.xq,B_Ki,Ki_cuant.s);
    fprintf("\n Los errores son:");
    fprintf("\n Error magnitud en 0: %d",Er_m_f0_real);
    fprintf("\n Error magnitud en fc: %d",Er_m_fc_real);
    fprintf("\n Error fase en fc: %d",Er_f_fc_real);
    break
end
end



figure('Name','TF de PID escalada y  cuantificado')

bodef(Gcz_cuant);

figure('Name','TF comparado')
hold on
bodef(Gcz_cuant);
bodef(Gcz_ncuant);

%Funcion B2C
function [wk] = Qn(x,n)
x1 = x;
neg = (x<0);
x = abs(x);
E = floor(log2(x));
F = 2^(log2(x)-E);
q = E-(n-2);
wd = round(F*2^(n-2));
if (wd==1)
wd = round(F*2^(n-3));
q = q+1;
end
if (neg)
xq = -wd*2^q;
s = [' 1',dec2bin(-wd+2^(n-1),n-1)];
wd = -wd;
else
xq = wd*2^q;
s = [' 0',dec2bin(wd,n-1)];
end;

wk.xq = xq;    %el valor qe se obtiene al final
wk.w = wd;      % Valor si no fuese con radix
wk.q = q;       %Escala
wk.n = n;       %Numero de bits
wk.s = s;        %valor en bits
wk.dx = xq-x1;  
return;
end

function h = bodef(x)
   P = bodeoptions; P.FreqUnits = 'Hz';
   h = bodeplot(x,P);
end

