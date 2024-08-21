clc;
clear ;
Var_Vo=5; %variacion en % maxima permitida en Vo, en SS.
Vo=1.8;
Vfs=2;       %Voltaje maximo que permite el ADC, el minimo se supone que es 0
Vref=1.8; %Que este en medio de los bins del ADC, no es necesario, pero NO debe estar en los extremos
H=1;
B_ADC=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
qvo_r=Vfs/(2^B_ADC);
fprintf('Valor de Vref cuantificado:')
disp(Qn(Vref,qvo_r,Vfs));
Vref_cuant_ini=Qn(Vref,qvo_r,Vfs);
B_min_ADC=ceil(log2( (Vfs*100)/(H*Vo*Var_Vo) ));
qvo=Vfs/(2^B_min_ADC);

while ((Qn(Vref,qvo,Vfs)+qvo) > (Vref+Var_Vo/100*Vref)*H) || ((Qn(Vref,qvo,Vfs)-qvo)<(Vref-Var_Vo/100*Vref)*H)
B_min_ADC=B_min_ADC+1;
qvo=Vfs/(2^B_min_ADC);
display(B_min_ADC)
display(qvo)
disp(Qn(Vref,qvo,Vfs))
end

Vref_cuant=Qn(Vref,qvo,Vfs);
display(B_min_ADC)
disp(Vref_cuant)
display(qvo)

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