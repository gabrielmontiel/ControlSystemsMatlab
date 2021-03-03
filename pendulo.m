M = .5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

states = {'x' 'x_p' 'phi' 'phi_p' };
inputs = {'u'};
outputs = {'x'; 'phi'};

sys_ss = ss(A,B,C,D,'statename',states,'inputname',inputs,'outputname',outputs)

%%Observabilidad y controlabilidad
%Controlabilidad
r1=[];
for i = 0:length(A)-1
    r1=horzcat(r1,A^(i)*B);
end
Controlabilidad=rank(r1)
%Observabilidad
r2=[];
for i = 0:length(A)-1
    r2=vertcat(r2,C*A^(i));
end
observabilidad = rank(r2)

%Estabilidad
e=eig(A)
%El sistema es estable si todos los eigenvalues son negativos
%El sistema no es estable ya que los polos son
%         0
%    -0.1428
%    -5.6041
%     5.5651 <--- este polo es positivo, inestable


sys_ss;
sys_tf=tf(sys_ss);
pzmap(sys_ss);


%% PID pole placement
%sys1 es sistema a lazo abierto.
%sys1_sym es el sistema en simbolico
%sysFb1 es el sistema con feedback PID simbolico
%pid_sym es el controlador PID en simbolico

sys1=sys_tf(2);

syms ki kd kp s
H=1; %Sensor
%s=tf('s'); recordatorio

[num,den]=tfdata(sys1,'v');

sys1_sym=poly2sym(num,s)/poly2sym(den,s);
pid_sym=kp+ki/s+kd*s;


sysPID_sym=sys1_sym*pid_sym;
sysFb1=sysPID_sym/(1+sysPID_sym*H);

[~,d]=numden(sysFb1);
d=collect(d,s);
coeficientes=fliplr(coeffs(d,s)); %Ordenando de mayor a menor
coeficientes=coeficientes/coeficientes(1); %Diviendo por la S mayor

poles=[-3,-8,-10]; % Ingresar POLOS
pFunction=1; %Iniciar
for i=poles
    pFunction=pFunction*(s-i); %Funcion con los Polos
end
pFunction=collect(pFunction);
coeffs_pid=fliplr(coeffs(pFunction)); %Ordenando de mayor a menor

%ya aqui estan Coeffs_pid y Coeficientes como vectores de s mayor a s menor
resultados = zeros(1,3);
for i = 2:length(coeffs_pid)
    x=coeficientes(i)==coeffs_pid(i);
    x=subs(x);
    y=isolate(x,symvar(x));
    value=solve(y,symvar(y));
    assignin('base',string(symvar(y)),value);
end
% 
% prueba=resultados;
% for i=1:length(coeffs_pid) %Calculando soluciones
%     if i == 1
%         continue %La primera no se hace es 1=1.
%     end
%     %pruebita(i)=isolate(coeficientes(i)==c(i),symvar(coeficientes(i)));
%     
%     resultados(i)=solve(coeficientes(i)==coeffs_pid(i));
%     assignin('base',string(symvar(coeficientes(i))),resultados(i));
% end

% Calculando Respuesta al escalon con el respectivo PID
ki=round(double(ki),3);
kp=round(double(kp),3);
kd=round(double(kd),3);

PID1=pid(kp,ki,kd);

syspid1=feedback(PID1*sys1,H);

%% prueba
[~,d23]=tfdata(syspid1,'v');
poly23=poly2sym(d23,s);

%% asignar polos PI
