clc
clear all

%solving for volumetric heating based on thermophysical properties of CPA
%heat capacity

T1 = -135.4;%deg.C
T2 = -80;
T3 = -40;
R1 = 1/60;%deg.C/min*1min/60*second
R2 = 100/60;%deg.C/second
dt = 12;%s

t1=(T2-T1)/R1;
t2=(T3-T2)/R2;
t = round(t1/dt)+round(t2/dt);
tm = linspace(1,t,t)';%column of timesteps

%column of target temperatures based on desired cooling rate
Texp = zeros(size(tm));
Rexp = zeros(size(tm));
Rexp(1)=R1;
Texp(1) = T1+(R1*dt);
for i=2:size(Texp)
if Texp(i-1)<=T2
    Texp(i)=Texp(i-1)+(R1*dt);
    Rexp(i)=R1;
else
    Texp(i)=Texp(i-1)+(R2*dt);
    Rexp(i)=R2;
end
end


cp =[-209.55	97.33076786
-194.95	272.0626701
-181.65	496.9166073
-171.25	696.8750556
-162.05	879.9495951
-150.25	1112.68269
-136.95	1361.637579
-120.85	1635.090689
-106.25	1853.919388
-90.75	2057.188406
-73.65	2249.201486
-64.25	2341.897455
-55.95	2417.378459
-43.95	2516.364512
-35.35	2580.589576
-30	2591.049982
-25	2622.267719
-20	2666.347092
-10	2737.976073
0	2798.585211
20	2925.313408];% J/kg-K

p = [-165	1125
20	1050];%kg/m3

%temp*cp=J/kg;
%J/kg*kg/m3=J/m3;
%J/m3 /(dt) =W/m3;
A=size(cp)
Cexp=zeros(size(Texp));
for i=1:size(Texp)
    Tr=Texp(i).*(ones(A(1),1))
    Td=abs(cp(:,1)-Tr);
    [M,I] = min(Td);
    Cexp(i) = cp(I,2);
end

cp =[-209.55	97.33076786
-194.95	272.0626701
-181.65	496.9166073
-171.25	696.8750556
-162.05	879.9495951
-150.25	1112.68269
-136.95	1361.637579
-120.85	1635.090689
-106.25	1853.919388
-90.75	2057.188406
-73.65	2249.201486
-64.25	2341.897455
-55.95	2417.378459
-43.95	2516.364512
-35.35	2580.589576
-30	2591.049982
-25	2622.267719
-20	2666.347092
-10	2737.976073
0	2798.585211
20	2925.313408];% J/kg-K

p = [-165	1125
20	1050];%kg/m3

pl = polyfit(p(:,1),p(:,2),1);

pexp = polyval(pl,Texp);

qexp = Cexp.*pexp.*Rexp;





