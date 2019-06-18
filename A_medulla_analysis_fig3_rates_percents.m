clc
%%
clear all
medl = 63205;
medm = 103574;
% corl = 13457;
% filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\A2_split_12p5_cor1'
% load(filename)
% cor = horzcat(cor{1},cor{2},cor{3},cor{4});
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\A2_split_12p5_med1'
load(filename)
medu = horzcat(med{1},med{2},med{3},med{4});

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\A2_split_12p5_med2'
load(filename)
med2 = horzcat(med{1},med{2},med{3},med{4});

% t = 12.5;
% filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonex.csv'
% load(filename)
% X = corsideonex(13458:48002);
% filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideoney.csv'
% load(filename)
% Y = corsideoney(13458:48002);
% filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonez.csv'
% load(filename)
% Z = corsideonez(13458:48002);
% 
% Ccoords = horzcat(X,Y,Z);

%%%medulla
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medonex.csv'
load(filename)
X = medonex(63206:65829);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medoney.csv'
load(filename)
Y = medoney(63206:65829);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medonez.csv'
load(filename)
Z = medonez(63206:65829);

Mcoords = horzcat(X,Y,Z);

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medtwox.csv'
load(filename)
X2 = medtwox(103575:106639);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medtwoy.csv'
load(filename)
Y2 = medtwoy(103575:106639);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medtwoz.csv'
load(filename)
Z2 = medtwoz(103575:106639);

Mcoords2 = horzcat(X2,Y2,Z2);

[q,r]=size(medu);
med1 = -22*ones(q,1);
medu = horzcat(med1,medu);
[q,r]=size(medu);
timestep = 12.5;%seconds
tminut = timestep/60;%minute



time = 12.5*linspace(1,880,880);
time = horzcat(0,time);
%points of interest

time = 12.5*linspace(1,880,880);
time = horzcat(0,time);

medR = ones(q,r-4);
medT = ones(q,r-4);

for i = 1:q    
    for j=3:r-2
        psb=horzcat(medu(i,j-2),medu(i,j-1),medu(i,j),medu(i,j+1),medu(i,j+2));
        tsb=time(j-2:j+2);
        sb = polyfit(tsb,psb,1);
        medR(i,j) = sb(1);    
    end
end

for i=1:q
    for j=3:r-2
        if medu(i,j)<=-40 && medu(i,j)>=-80;
            medT(i,j) = 0;
        end
    end
end

Z = ones(q,1);
medT = horzcat(Z,Z,medT);
medP = medT.*medR;
medY = ones(q,1);
[s,t]=size(medP);

for i=1:s
    for j=1:t
    if medP(i,j)==0
        medY(i)=0;
    end
    end
end

[q,r]=size(med2);
med12 = -22*ones(q,1);
med2 = horzcat(med12,med2);
[q,r]=size(med2);
timestep = 12.5;%seconds
tminut = timestep/60;%minute



time = 12.5*linspace(1,880,880);
time = horzcat(0,time);
%points of interest

time = 12.5*linspace(1,880,880);
time = horzcat(0,time);

medR2 = ones(q,r-4);
medT2 = ones(q,r-4);

for i = 1:q    
    for j=3:r-2
        psb=horzcat(med2(i,j-2),med2(i,j-1),med2(i,j),med2(i,j+1),med2(i,j+2));
        tsb=time(j-2:j+2);
        sb = polyfit(tsb,psb,1);
        medR2(i,j) = sb(1);    
    end
end

for i=1:q
    for j=3:r-2
        if med2(i,j)<=-40 && med2(i,j)>=-80;
            medT2(i,j) = 0;
        end
    end
end

Z2 = ones(q,1);
medT2 = horzcat(Z2,Z2,medT2);
medP2 = medT2.*medR2;
medY2 = ones(q,1);
[s,t]=size(medP2);

for i=1:s
    for j=1:t
    if medP2(i,j)==0
        medY2(i)=0;
    end
    end
end

%chunks instead of individual lines
%
%corR=-corR*60;
medR=-medR*60;
figure(3)
hold on
plot(cor(b,3:end-2),corR(1,3:end),'-c','LineWidth',2)
hold on
plot(cor(c,3:end-2),corR(2,3:end),'-c','LineWidth',2)
plot(cor(d,3:end-2),corR(3,3:end),'-c','LineWidth',2)
plot(cor(e,3:end-2),corR(4,3:end),'-c','LineWidth',2)
plot(cor(f,3:end-2),corR(5,3:end),'-c','LineWidth',2)
plot(cor(k,3:end-2),corR(6,3:end),'-c','LineWidth',2)
plot(cor(f,3:end-2),(0.1*ones(size(corR(5,3:end)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,3:end-2),medR(1,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(time(3:end-2),medR(2,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(g,3:end-2),medR(2,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(h,3:end-2),medR(3,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(m,3:end-2),medR(4,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(n,3:end-2),medR(5,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(p,3:end-2),medR(6,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(cor(f,3:end-2),(7*ones(size(corR(5,3:end)))),'--k','LineWidth',2)
hold off

