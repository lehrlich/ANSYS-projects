clear all
clc

%%
clc
clear all

medl = 63205;
corl = 13457;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\B2_split_12p5_cor1'
load(filename)
cor = horzcat(cor{1},cor{2},cor{3},cor{4});
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\B2_split_12p5_med1'
load(filename)
med = horzcat(med{1},med{2},med{3},med{4});
t = 12.5;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonex.csv'
load(filename)
X = corsideonex(13458:48002);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideoney.csv'
load(filename)
Y = corsideoney(13458:48002);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonez.csv'
load(filename)
Z = corsideonez(13458:48002);

Ccoords = horzcat(X,Y,Z);
sliceC = zeros(size(Ccoords));

for i=1:length(Ccoords)
    if Ccoords(i,1)>=-5.334e-3 & Ccoords(i,1)<=-5.0694e-3
        sliceC(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceC(i,:) = NaN;
    end
end
ccoord_slice = Ccoords.*sliceC;

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
sliceM = zeros(size(Mcoords));

for i=1:length(Mcoords)
    if Mcoords(i,1)>=-5.334e-3 & Mcoords(i,1)<=-5.0694e-3
        sliceM(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceM(i,:) = NaN;
    end
end
mcoord_slice = Mcoords.*sliceM;

%calculate rates using central differences, etc.
[m,n]=size(cor);
cor1 = -22*ones(m,1);
cor = horzcat(cor1,cor);
[m,n]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute
for j=1:n-1
    curt = cor(:,j);
    next = cor(:,j+1);
    corR(:,j) = ((curt-next)/tminut);
end

[m,n]=size(med);
med1 = -22*ones(m,1);
med = horzcat(med1,med);
[m,n]=size(med);

for j=1:n-1
    curt = med(:,j);
    next = med(:,j+1);
    medR(:,j) = ((curt-next)/tminut);
end


    
    
%points of interest
a = 65422;
b = 13502;
c = 13515;
d = 13485;
e = 13536;
f = 14061;
g = 64746;
h = 63358;
k = 13467;
m = 63240;
n = 63212;
p = 63261;

%cor points
b = b - corl;
c = c - corl;
d = d - corl;
e = e - corl;
f = f - corl;
k = k - corl;

%med points
a = a - medl;
g = g - medl;
h = h - medl;
m = m - medl;
n = n - medl;
p = p - medl;

time = 12.5*linspace(1,400,400);
time = horzcat(0,time);


figure(3)
plot(cor(b,2:end),corR(b,:),'-b','LineWidth',2)
hold on
plot(cor(c,2:end),corR(c,:),'-b','LineWidth',2)
plot(cor(d,2:end),corR(d,:),'-b','LineWidth',2)
plot(cor(e,2:end),corR(e,:),'-b','LineWidth',2)
plot(cor(f,2:end),corR(f,:),'-b','LineWidth',2)
plot(cor(k,2:end),corR(k,:),'-b','LineWidth',2)
plot(cor(f,2:end),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,2:end),medR(a,:),'-r','LineWidth',2)
plot(med(g,2:end),medR(g,:),'-r','LineWidth',2)
plot(med(h,2:end),medR(h,:),'-r','LineWidth',2)
plot(med(m,2:end),medR(m,:),'-r','LineWidth',2)
plot(med(n,2:end),medR(n,:),'-r','LineWidth',2)
plot(med(p,2:end),medR(p,:),'-r','LineWidth',2)
plot(cor(f,2:end),(7*ones(size(corR(f,:)))),'--k','LineWidth',2)
hold off




%%
clear all
clc

medl = 31441;
corl = 5847;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\rA2_split_12p5_cor1'
load(filename)
cor = horzcat(cor{1},cor{2},cor{3},cor{4});
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\rA2_split_12p5_med1'
load(filename)
med = horzcat(med{1},med{2},med{3},med{4});
t = 12.5;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonex.csv'
load(filename)
X = corsideonex(5848:22425);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideoney.csv'
load(filename)
Y = corsideoney(5848:22425);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonez.csv'
load(filename)
Z = corsideonez(5848:22425);

Ccoords = horzcat(X,Y,Z);
sliceC = zeros(size(Ccoords));

for i=1:length(Ccoords)
     if Ccoords(i,1)>=-1.9203e-3 & Ccoords(i,1)<=-1.7e-3
        sliceC(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceC(i,:) = NaN;
    end
end
ccoord_slice = Ccoords.*sliceC;

%%%medulla
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medonex.csv'
load(filename)
X = medonex(31442:33115);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medoney.csv'
load(filename)
Y = medoney(31442:33115);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medonez.csv'
load(filename)
Z = medonez(31442:33115);

Mcoords = horzcat(X,Y,Z);
sliceM = zeros(size(Mcoords));

for i=1:length(Mcoords)
    if Ccoords(i,1)>=-1.9203e-3 & Ccoords(i,1)<=-1.7e-3
        sliceM(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceM(i,:) = NaN;
    end
end
mcoord_slice = Mcoords.*sliceM;
%calculate rates using central differences, etc.
[m,n]=size(cor);
cor1 = -22*ones(m,1);
cor = horzcat(cor1,cor);
[m,n]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute
for j=1:n-1
    curt = cor(:,j);
    next = cor(:,j+1);
    corR(:,j) = ((curt-next)/tminut);
end

[m,n]=size(med);
med1 = -22*ones(m,1);
med = horzcat(med1,med);
[m,n]=size(med);

for j=1:n-1
    curt = med(:,j);
    next = med(:,j+1);
    medR(:,j) = ((curt-next)/tminut);
end    
    
%points of interest
a = 32887;
b = 5889;
c = 5901;
d = 5872;
e = 5922;
f = 6236;
g = 31549;
h = 31547;
k = 5856;
m = 31912;
n = 31448;
p = 31481;

%cor points
b = b - corl;
c = c - corl;
d = d - corl;
e = e - corl;
f = f - corl;
k = k - corl;

%med points
a = a - medl;
g = g - medl;
h = h - medl;
m = m - medl;
n = n - medl;
p = p - medl;

time = 12.5*linspace(1,300,300);
time = horzcat(0,time);


%chunks instead of individual lines
%cortex
figure(3)
hold on
plot(cor(b,2:end),corR(b,:),'-g','LineWidth',2)
hold on
plot(cor(c,2:end),corR(c,:),'-g','LineWidth',2)
plot(cor(d,2:end),corR(d,:),'-g','LineWidth',2)
plot(cor(e,2:end),corR(e,:),'-g','LineWidth',2)
plot(cor(f,2:end),corR(f,:),'-g','LineWidth',2)
plot(cor(k,2:end),corR(k,:),'-g','LineWidth',2)
plot(cor(f,2:end),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)

figure(3)
hold on
plot(med(a,2:end),medR(a,:),'-m','LineWidth',2)
plot(med(g,2:end),medR(g,:),'-m','LineWidth',2)
plot(med(h,2:end),medR(h,:),'-m','LineWidth',2)
plot(med(m,2:end),medR(m,:),'-m','LineWidth',2)
plot(med(n,2:end),medR(n,:),'-m','LineWidth',2)
plot(med(p,2:end),medR(p,:),'-m','LineWidth',2)
plot(cor(f,2:end),(7*ones(size(corR(f,:)))),'--k','LineWidth',2)

%%
clear all
%%
clear all
clc

medl = 31441;
corl = 5847;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\rB2_split_12p5_cor1'
load(filename)
cor = horzcat(cor{1},cor{2},cor{3},cor{4});
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\rB2_split_12p5_med1'
load(filename)
med = horzcat(med{1},med{2},med{3},med{4});
t = 12.5;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonex.csv'
load(filename)
X = corsideonex(5848:22425);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideoney.csv'
load(filename)
Y = corsideoney(5848:22425);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonez.csv'
load(filename)
Z = corsideonez(5848:22425);

Ccoords = horzcat(X,Y,Z);
sliceC = zeros(size(Ccoords));

for i=1:length(Ccoords)
     if Ccoords(i,1)>=-1.9203e-3 & Ccoords(i,1)<=-1.7e-3
        sliceC(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceC(i,:) = NaN;
    end
end
ccoord_slice = Ccoords.*sliceC;

%%%medulla
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medonex.csv'
load(filename)
X = medonex(31442:33115);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medoney.csv'
load(filename)
Y = medoney(31442:33115);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\medonez.csv'
load(filename)
Z = medonez(31442:33115);

Mcoords = horzcat(X,Y,Z);
sliceM = zeros(size(Mcoords));

for i=1:length(Mcoords)
    if Ccoords(i,1)>=-1.9203e-3 & Ccoords(i,1)<=-1.7e-3
        sliceM(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceM(i,:) = NaN;
    end
end
mcoord_slice = Mcoords.*sliceM;
%calculate rates using central differences, etc.
[m,n]=size(cor);
cor1 = -22*ones(m,1);
cor = horzcat(cor1,cor);
[m,n]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute
for j=1:n-1
    curt = cor(:,j);
    next = cor(:,j+1);
    corR(:,j) = ((curt-next)/tminut);
end

[m,n]=size(med);
med1 = -22*ones(m,1);
med = horzcat(med1,med);
[m,n]=size(med);

for j=1:n-1
    curt = med(:,j);
    next = med(:,j+1);
    medR(:,j) = ((curt-next)/tminut);
end    
    
%points of interest
a = 32887;
b = 5889;
c = 5901;
d = 5872;
e = 5922;
f = 6236;
g = 31549;
h = 31547;
k = 5856;
m = 31912;
n = 31448;
p = 31481;

%cor points
b = b - corl;
c = c - corl;
d = d - corl;
e = e - corl;
f = f - corl;
k = k - corl;

%med points
a = a - medl;
g = g - medl;
h = h - medl;
m = m - medl;
n = n - medl;
p = p - medl;

time = 12.5*linspace(1,300,300);
time = horzcat(0,time);

figure(3)
hold on
plot(cor(b,2:end),corR(b,:),'-k','LineWidth',2)
hold on
plot(cor(c,2:end),corR(c,:),'-k','LineWidth',2)
plot(cor(d,2:end),corR(d,:),'-k','LineWidth',2)
plot(cor(e,2:end),corR(e,:),'-k','LineWidth',2)
plot(cor(f,2:end),corR(f,:),'-k','LineWidth',2)
plot(cor(k,2:end),corR(k,:),'-k','LineWidth',2)
plot(cor(f,2:end),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)

figure(3)
hold on
plot(med(a,2:end),medR(a,:),'-y','LineWidth',2)
plot(med(g,2:end),medR(g,:),'-y','LineWidth',2)
plot(med(h,2:end),medR(h,:),'-y','LineWidth',2)
plot(med(m,2:end),medR(m,:),'-y','LineWidth',2)
plot(med(n,2:end),medR(n,:),'-y','LineWidth',2)
plot(med(p,2:end),medR(p,:),'-y','LineWidth',2)
plot(cor(f,2:end),(7*ones(size(corR(f,:)))),'--k','LineWidth',2)

%%
clear all
medl = 63205;
corl = 13457;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\A2_split_12p5_cor1'
load(filename)
cor = horzcat(cor{1},cor{2},cor{3},cor{4});
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\Compiled_split_results_12p5\A2_split_12p5_med1'
load(filename)
med = horzcat(med{1},med{2},med{3},med{4});
t = 12.5;
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonex.csv'
load(filename)
X = corsideonex(13458:48002);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideoney.csv'
load(filename)
Y = corsideoney(13458:48002);
filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\corsideonez.csv'
load(filename)
Z = corsideonez(13458:48002);

Ccoords = horzcat(X,Y,Z);
sliceC = zeros(size(Ccoords));

for i=1:length(Ccoords)
    if Ccoords(i,1)>=-5.334e-3 & Ccoords(i,1)<=-5.0694e-3
        sliceC(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceC(i,:) = NaN;
    end
end
ccoord_slice = Ccoords.*sliceC;

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
sliceM = zeros(size(Mcoords));

for i=1:length(Mcoords)
    if Mcoords(i,1)>=-5.334e-3 & Mcoords(i,1)<=-5.0694e-3
        sliceM(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceM(i,:) = NaN;
    end
end
mcoord_slice = Mcoords.*sliceM;

[m,n]=size(cor);
cor1 = -22*ones(m,1);
cor = horzcat(cor1,cor);
[m,n]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute
for j=1:n-1
    curt = cor(:,j);
    next = cor(:,j+1);
    corR(:,j) = ((curt-next)/tminut);
end

[m,n]=size(med);
med1 = -22*ones(m,1);
med = horzcat(med1,med);
[m,n]=size(med);

for j=1:n-1
    curt = med(:,j);
    next = med(:,j+1);
    medR(:,j) = ((curt-next)/tminut);
end


time = 12.5*linspace(1,880,880);
time = horzcat(0,time);
%points of interest
a = 65422;
b = 13502;
c = 13515;
d = 13485;
e = 13536;
f = 14061;
g = 64746;
h = 63358;
k = 13467;
m = 63240;
n = 63212;
p = 63261;

%cor points
b = b - corl;
c = c - corl;
d = d - corl;
e = e - corl;
f = f - corl;
k = k - corl;

%med points
a = a - medl;
g = g - medl;
h = h - medl;
m = m - medl;
n = n - medl;
p = p - medl;

time = 12.5*linspace(1,880,880);
time = horzcat(0,time);
%%


%chunks instead of individual lines
%

figure(3)
hold on
plot(cor(b,2:end),corR(b,:),'-c','LineWidth',2)
hold on
plot(cor(c,2:end),corR(c,:),'-c','LineWidth',2)
plot(cor(d,2:end),corR(d,:),'-c','LineWidth',2)
plot(cor(e,2:end),corR(e,:),'-c','LineWidth',2)
plot(cor(f,2:end),corR(f,:),'-c','LineWidth',2)
plot(cor(k,2:end),corR(k,:),'-c','LineWidth',2)
plot(cor(f,2:end),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)

figure(3)
hold on
plot(med(a,2:end),medR(a,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(g,2:end),medR(g,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(h,2:end),medR(h,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(m,2:end),medR(m,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(n,2:end),medR(n,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(p,2:end),medR(p,:),'-','color',[1 0.5 0],'LineWidth',2)
plot(cor(f,2:end),(7*ones(size(corR(f,:)))),'--k','LineWidth',2)
