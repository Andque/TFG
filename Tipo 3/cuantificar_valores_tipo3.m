clc
close all
clear all
%hay que ver en fc y en casi 0Hz la diferencia en % y dependiendo de el
%maximo permitido se escoge un numero de bits u otro.

%Valores constantes
V_max_ADC=1;
V_min_ADC=0;
B_ADC=12;
Nr=2047; %Poner el NR de verdad que en este paso se usa 
K1=0.014021043256821 ;
K2=59.688132747676260;
K3=-22.252867694489982;
K4=-30.083066311429160;
C1=1;
C2=0.487004818221354;
C3=0.487004810535232;

Ts=2.5e-05; %Tiempo Tsw
fc= (1/Ts) * (1/60); %frecuencia de corte

%%Especifico de este paso

B_min_r=4  ;  %minimo numero de bits para los valores del regulador MINIMO 2!!!
B_max_r=10; %maximo numero de bits permitidos
Er_m_f0=5; %Porcentaje de diferencia en la magnitud en frecuancia 0
Er_f_fc=1;  %Porcentaje de diferencia en la fase en frecuancia fc
Er_m_fc=0.56;  %Porcentaje de diferencia en la magnitud en frecuancia fc
boolean Ok=false;
%Tipo 3 en // con el escalado!! Que en este caso es multiplicar por Nr y qvs
qvs=(V_max_ADC-V_min_ADC)/(2^B_ADC);
K1_e=K1*Nr*qvs;
K2_e=K2*Nr*qvs;
K3_e=K3*Nr*qvs;
K4_e=K4*Nr*qvs;
C1_e=C1;
C2_e=C2;
C32_e=-(C3^2); %El valor de C3^2
C31_e=C3*2; %El valor de C3*2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Calculo de bits minimos   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Introducir la funcion de trasnferencia del PID sin cuantificar y con el
%escalado!!!!
z = tf('z',Ts);
Gcz3_pe=(K1_e/(1-C1_e*z^-1))+(K2_e/(1-C2_e*z^-1))+K4_e + (K3_e/(-C32_e*(z^-2)-C31_e*(z^-1)+1));
figure('Name','TF de tipo 3 escalado y sin cuantificar')
hold on
bodef(Gcz3_pe);
[mag_ncuant_0,phase_ncuant_0] = bode(Gcz3_pe,0.1);
[mag_ncuant_fc,phase_ncuant_fc] = bode(Gcz3_pe,fc*2*pi);

%calcular con los distintos valores de Kp,Ki y Kd cuantificados cual es el
%k hace que el error sea el maximo permitido con el minimo numero de bits.

for B_K1 = B_min_r:B_max_r

K1_cuant=Qn(K1_e,B_K1); %Cuantificar K1 con los bits correspondientes
        for B_K2 = B_min_r:B_max_r
         K2_cuant=Qn(K2_e,B_K2);%Cuantificar K2 con los bits correspondientes
            for B_K3= B_min_r:B_max_r
            K3_cuant=Qn(K3_e,B_K3);%Cuantificar K con los bits correspondientes    
                for B_K4=B_min_r:B_max_r
                K4_cuant=Qn(K4_e,B_K4);%Cuantificar K con los bits correspondientes
                    for B_C1=3:3
                    C1_cuant=Qn(C1_e,B_C1);%Cuantificar C con los bits correspondientes
                        for B_C2=B_min_r:B_max_r
                        C2_cuant=Qn(C2_e,B_C2);%Cuantificar C con los bits correspondientes    
                            for B_C31=B_min_r:B_max_r
                            C31_cuant=Qn(C31_e,B_C31);%Cuantificar C con los bits correspondientes
                                for B_C32=B_min_r:B_max_r
                                C32_cuant=Qn(C32_e,B_C32);%Cuantificar C con los bits correspondientes
                                % Crear una nueva TF con los nuevos valores
                                z = tf('z',Ts);
                                Gcz3_cuant=(K1_cuant.xq/(1-C1_cuant.xq*z^-1))+(K2_cuant.xq/(1-C2_cuant.xq*z^-1))+K4_cuant.xq + (K3_cuant.xq/(-C32_cuant.xq*(z^-2)-C31_cuant.xq*(z^-1)+1));

                                [mag_cuant_0,phase_cuant_0] = bode(Gcz3_cuant,0.1);
                                [mag_cuant_fc,phase_cuant_fc] = bode(Gcz3_cuant,fc*2*pi);

                                % Ver si se cumplen las 3 condiciones
                                Er_m_f0_real=abs (abs(mag_ncuant_0-mag_cuant_0) / mag_ncuant_0) *100;
                                Er_m_fc_real=abs( abs(mag_ncuant_fc-mag_cuant_fc) / mag_ncuant_fc * 100);
                                Er_f_fc_real=abs( abs(phase_ncuant_fc-phase_cuant_fc) / phase_ncuant_fc * 100);
                                fprintf("\nB_K1=%d, B_K2=%d, B_K3=%d, B_K4=%d, B_C1=%d, B_C2=%d, B_C31=%d, B_C32=%d",B_K1,B_K2,B_K3,B_K4,B_C1,B_C2,B_C31,B_C32);
                                Ok=false;
                                if( Er_m_f0_real<Er_m_f0 && Er_m_fc_real<Er_m_fc && Er_f_fc_real<Er_f_fc  )
                                Ok=true;
                                break
                                end 
                                end
                            if( Ok )
                            break
                            end 
                            end
                            
                        if( Ok)
                        break
                        end     
                        end

                    if( Ok)
                    break
                    end     
                    end

                if( Ok )
                break
                end 
                end

            if( Ok  )
            break
            end     
            end

        if( Ok  )
        break
        end 
        end

if( Ok  )
break
end 
end

         





    if( Ok )
        fprintf("\n El valor para que se cumpla los errores maximos permitidos son");
        fprintf("\n K1=  %d B_K1= %d",K1_cuant.xq,B_K1);
        fprintf("\n K2=  %d B_K2= %d",K2_cuant.xq,B_K2);
        fprintf("\n K3=  %d B_K3= %d",K3_cuant.xq,B_K3);
        fprintf("\n K4=  %d B_K4= %d",K4_cuant.xq,B_K4);
        fprintf("\n C1=  %d B_C1= %d",C1_cuant.xq,B_C1);
        fprintf("\n C2=  %d B_C2= %d",C2_cuant.xq,B_C2);
        fprintf("\n C31=  %d B_C31= %d",C31_cuant.xq,B_C31);
        fprintf("\n C32=  %d B_C32= %d",C32_cuant.xq,B_C32);

        fprintf("\n Los errores son:");
        fprintf("\n Error magnitud en 0: %d",Er_m_f0_real);
        fprintf("\n Error magnitud en fc: %d",Er_m_fc_real);
        fprintf("\n Error fase en fc: %d",Er_f_fc_real);
        
    end
  


figure('Name','TF de PID escalada y  cuantificado')

bodef(Gcz3_cuant);

figure('Name','TF comparado')
hold on
bodef(Gcz3_cuant);
bodef(Gcz3_pe);

%Funcion B2C
function [wk] = Qn(x,n)
x1 = x;
neg= (x<0);
x = abs(x);
E = floor(log2(x));
F = 2^(log2(x)-E);
q = E-(n-2);
wd = round(F*2^(n-2));
if (wd==1)
wd = round(F*2^(n-3));
q = q+1;
end;
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