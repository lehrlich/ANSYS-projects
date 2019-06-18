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


%calculate rates using central differences, etc.
[q,r]=size(cor);
cor1 = -22*ones(q,1);
cor = horzcat(cor1,cor);
[q,r]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute

time = 12.5*linspace(1,400,400);
time = horzcat(0,time);

for j=3:r-2
    
    psb=horzcat(cor(b,j-2),cor(b,j-1),cor(b,j),cor(b,j+1),cor(b,j+2));
    tsb=time(j-2:j+2);
    sb = polyfit(tsb,psb,1);
    corR(1,j) = sb(1);
    
    psc=horzcat(cor(c,j-2),cor(c,j-1),cor(c,j),cor(c,j+1),cor(c,j+2));
    sc = polyfit(tsb,psc,1);
    corR(2,j) = sc(1);
    
    psd=horzcat(cor(d,j-2),cor(d,j-1),cor(d,j),cor(d,j+1),cor(d,j+2));
    sd = polyfit(tsb,psd,1);
    corR(3,j) = sd(1);
    
    pse=horzcat(cor(e,j-2),cor(e,j-1),cor(e,j),cor(e,j+1),cor(e,j+2));
    se = polyfit(tsb,pse,1);
    corR(4,j) = se(1);
    
    psf=horzcat(cor(f,j-2),cor(f,j-1),cor(f,j),cor(f,j+1),cor(f,j+2));
    sf = polyfit(tsb,psf,1);
    corR(5,j) = sf(1);
    
    psk= horzcat(cor(k,j-2),cor(k,j-1),cor(k,j),cor(k,j+1),cor(k,j+2));
    sk = polyfit(tsb,psk,1);
    corR(6,j) = sk(1);
    
end

[q,r]=size(med);
med1 = -22*ones(q,1);
med = horzcat(med1,med);
[q,r]=size(med);

for j=3:r-2
    psa=horzcat(med(a,j-2),med(a,j-1),med(a,j),med(a,j+1),med(a,j+2));
    tsb=time(j-2:j+2);
    sa = polyfit(tsb,psa,1);
    medR(1,j) = sa(1);
    
    psg=horzcat(med(g,j-2),med(g,j-1),med(g,j),med(g,j+1),med(g,j+2));
    sg = polyfit(tsb,psg,1);
    medR(2,j) = sg(1);
    
    psh=horzcat(med(h,j-2),med(h,j-1),med(h,j),med(h,j+1),med(h,j+2));
    sh = polyfit(tsb,psh,1);
    medR(3,j) = sh(1);
    
    psm=horzcat(med(m,j-2),med(m,j-1),med(m,j),med(m,j+1),med(m,j+2));
    sm = polyfit(tsb,psm,1);
    medR(4,j) = sm(1);
    
    psn=horzcat(med(n,j-2),med(n,j-1),med(n,j),med(n,j+1),med(n,j+2));
    sn = polyfit(tsb,psn,1);
    medR(5,j) = sn(1);
    
    psp=horzcat(med(p,j-2),med(p,j-1),med(p,j),med(p,j+1),med(p,j+2));
    sp = polyfit(tsb,psp,1);
    medR(6,j) = sp(1);
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

corR=-corR*60;
medR=-medR*60;
figure(3)
plot(cor(b,3:end-2),corR(1,3:end),'-b','LineWidth',2)
hold on
plot(cor(c,3:end-2),corR(2,3:end),'-b','LineWidth',2)
plot(cor(d,3:end-2),corR(3,3:end),'-b','LineWidth',2)
plot(cor(e,3:end-2),corR(4,3:end),'-b','LineWidth',2)
plot(cor(f,3:end-2),corR(5,3:end),'-b','LineWidth',2)
plot(cor(k,3:end-2),corR(6,3:end),'-b','LineWidth',2)
plot(cor(f,3:end-2),(0.1*ones(size(corR(5,3:end)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,3:end-2),medR(1,3:end),'-r','LineWidth',2)
plot(med(g,3:end-2),medR(2,3:end),'-r','LineWidth',2)
plot(med(h,3:end-2),medR(3,3:end),'-r','LineWidth',2)
plot(med(m,3:end-2),medR(4,3:end),'-r','LineWidth',2)
plot(med(n,3:end-2),medR(5,3:end),'-r','LineWidth',2)
plot(med(p,3:end-2),medR(6,3:end),'-r','LineWidth',2)
plot(cor(f,3:end-2),(7*ones(size(corR(5,3:end)))),'--k','LineWidth',2)
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
[q,r]=size(cor);
cor1 = -22*ones(q,1);
cor = horzcat(cor1,cor);
[q,r]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute

    
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

for j=1:11
     corr=(cor(b,j+1)-cor(b,j))/12.5;
     corR(1,j) = corr;

     corr=(cor(c,j+1)-cor(c,j))/12.5;
     corR(2,j) = corr;
     
     corr=(cor(d,j+1)-cor(d,j))/12.5;
     corR(3,j) = corr;
     
     corr=(cor(e,j+1)-cor(e,j))/12.5;
     corR(4,j) = corr;

     corr=(cor(f,j+1)-cor(f,j))/12.5;
     corR(5,j) = corr;
     
     corr=(cor(k,j+1)-cor(k,j))/12.5;
     corR(6,j) = corr;
end
for j=11:r-2
    
    psb=horzcat(cor(b,j-2),cor(b,j-1),cor(b,j),cor(b,j+1),cor(b,j+2));
    tsb=time(j-2:j+2);
    sb = polyfit(tsb,psb,1);
    corR(1,j) = sb(1);
    
    psc=horzcat(cor(c,j-2),cor(c,j-1),cor(c,j),cor(c,j+1),cor(c,j+2));
    sc = polyfit(tsb,psc,1);
    corR(2,j) = sc(1);
    
    psd=horzcat(cor(d,j-2),cor(d,j-1),cor(d,j),cor(d,j+1),cor(d,j+2));
    sd = polyfit(tsb,psd,1);
    corR(3,j) = sd(1);
    
    pse=horzcat(cor(e,j-2),cor(e,j-1),cor(e,j),cor(e,j+1),cor(e,j+2));
    se = polyfit(tsb,pse,1);
    corR(4,j) = se(1);
    
    psf=horzcat(cor(f,j-2),cor(f,j-1),cor(f,j),cor(f,j+1),cor(f,j+2));
    sf = polyfit(tsb,psf,1);
    corR(5,j) = sf(1);
    
    psk= horzcat(cor(k,j-2),cor(k,j-1),cor(k,j),cor(k,j+1),cor(k,j+2));
    sk = polyfit(tsb,psk,1);
    corR(6,j) = sk(1);
    
end

[q,r]=size(med);
med1 = -22*ones(q,1);
med = horzcat(med1,med);
[q,r]=size(med);

for j=3:r-2
    psa=horzcat(med(a,j-2),med(a,j-1),med(a,j),med(a,j+1),med(a,j+2));
    sa = polyfit(tsb,psa,1);
    medR(1,j) = sa(1);
    
    psg=horzcat(med(g,j-2),med(g,j-1),med(g,j),med(g,j+1),med(g,j+2));
    sg = polyfit(tsb,psg,1);
    medR(2,j) = sg(1);
    
    psh=horzcat(med(h,j-2),med(h,j-1),med(h,j),med(h,j+1),med(h,j+2));
    sh = polyfit(tsb,psh,1);
    medR(3,j) = sh(1);
    
    psm=horzcat(med(m,j-2),med(m,j-1),med(m,j),med(m,j+1),med(m,j+2));
    sm = polyfit(tsb,psm,1);
    medR(4,j) = sm(1);
    
    psn=horzcat(med(n,j-2),med(n,j-1),med(n,j),med(n,j+1),med(n,j+2));
    sn = polyfit(tsb,psn,1);
    medR(5,j) = sn(1);
    
    psp=horzcat(med(p,j-2),med(p,j-1),med(p,j),med(p,j+1),med(p,j+2));
    sp = polyfit(tsb,psp,1);
    medR(6,j) = sp(1);
end

%chunks instead of individual lines
%cortex
corR=-corR*60;
medR=-medR*60;
figure(3)
hold on
plot(cor(b,1:end-2),corR(1,1:end),'-g','LineWidth',2)
hold on
plot(cor(c,1:end-2),corR(2,1:end),'-g','LineWidth',2)
plot(cor(d,1:end-2),corR(3,1:end),'-g','LineWidth',2)
plot(cor(e,1:end-2),corR(4,1:end),'-g','LineWidth',2)
plot(cor(f,1:end-2),corR(5,1:end),'-g','LineWidth',2)
plot(cor(k,1:end-2),corR(6,1:end),'-g','LineWidth',2)
%plot(cor(f,1:end-2),(0.1*ones(size(corR(5,3:end)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,3:end-2),medR(1,3:end),'-m','LineWidth',2)
plot(med(g,3:end-2),medR(2,3:end),'-m','LineWidth',2)
plot(med(h,3:end-2),medR(3,3:end),'-m','LineWidth',2)
plot(med(m,3:end-2),medR(4,3:end),'-m','LineWidth',2)
plot(med(n,3:end-2),medR(5,3:end),'-m','LineWidth',2)
plot(med(p,3:end-2),medR(6,3:end),'-m','LineWidth',2)
plot(cor(f,3:end-2),(7*ones(size(corR(5,3:end)))),'--k','LineWidth',2)
hold off



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
[q,r]=size(cor);
cor1 = -22*ones(q,1);
cor = horzcat(cor1,cor);
[q,r]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute

    
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

for j=3:r-2
    
    psb=horzcat(cor(b,j-2),cor(b,j-1),cor(b,j),cor(b,j+1),cor(b,j+2));
    tsb=time(j-2:j+2);
    sb = polyfit(tsb,psb,1);
    corR(1,j) = sb(1);
    
    psc=horzcat(cor(c,j-2),cor(c,j-1),cor(c,j),cor(c,j+1),cor(c,j+2));
    sc = polyfit(tsb,psc,1);
    corR(2,j) = sc(1);
    
    psd=horzcat(cor(d,j-2),cor(d,j-1),cor(d,j),cor(d,j+1),cor(d,j+2));
    sd = polyfit(tsb,psd,1);
    corR(3,j) = sd(1);
    
    pse=horzcat(cor(e,j-2),cor(e,j-1),cor(e,j),cor(e,j+1),cor(e,j+2));
    se = polyfit(tsb,pse,1);
    corR(4,j) = se(1);
    
    psf=horzcat(cor(f,j-2),cor(f,j-1),cor(f,j),cor(f,j+1),cor(f,j+2));
    sf = polyfit(tsb,psf,1);
    corR(5,j) = sf(1);
    
    psk= horzcat(cor(k,j-2),cor(k,j-1),cor(k,j),cor(k,j+1),cor(k,j+2));
    sk = polyfit(tsb,psk,1);
    corR(6,j) = sk(1);
    
end

[q,r]=size(med);
med1 = -22*ones(q,1);
med = horzcat(med1,med);
[q,r]=size(med);

for j=3:r-2
    psa=horzcat(med(a,j-2),med(a,j-1),med(a,j),med(a,j+1),med(a,j+2));
    sa = polyfit(tsb,psa,1);
    medR(1,j) = sa(1);
    
    psg=horzcat(med(g,j-2),med(g,j-1),med(g,j),med(g,j+1),med(g,j+2));
    sg = polyfit(tsb,psg,1);
    medR(2,j) = sg(1);
    
    psh=horzcat(med(h,j-2),med(h,j-1),med(h,j),med(h,j+1),med(h,j+2));
    sh = polyfit(tsb,psh,1);
    medR(3,j) = sh(1);
    
    psm=horzcat(med(m,j-2),med(m,j-1),med(m,j),med(m,j+1),med(m,j+2));
    sm = polyfit(tsb,psm,1);
    medR(4,j) = sm(1);
    
    psn=horzcat(med(n,j-2),med(n,j-1),med(n,j),med(n,j+1),med(n,j+2));
    sn = polyfit(tsb,psn,1);
    medR(5,j) = sn(1);
    
    psp=horzcat(med(p,j-2),med(p,j-1),med(p,j),med(p,j+1),med(p,j+2));
    sp = polyfit(tsb,psp,1);
    medR(6,j) = sp(1);
end

corR=-corR*60;
medR=-medR*60;
figure(3)
hold on
plot(cor(b,3:end-2),corR(1,3:end),'-k','LineWidth',2)
hold on
plot(cor(c,3:end-2),corR(2,3:end),'-k','LineWidth',2)
plot(cor(d,3:end-2),corR(3,3:end),'-k','LineWidth',2)
plot(cor(e,3:end-2),corR(4,3:end),'-k','LineWidth',2)
plot(cor(f,3:end-2),corR(5,3:end),'-k','LineWidth',2)
plot(cor(k,3:end-2),corR(6,3:end),'-k','LineWidth',2)
plot(cor(f,3:end-2),(0.1*ones(size(corR(5,3:end)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,3:end-2),medR(1,3:end),'-y','LineWidth',2)
plot(med(g,3:end-2),medR(2,3:end),'-y','LineWidth',2)
plot(med(h,3:end-2),medR(3,3:end),'-y','LineWidth',2)
plot(med(m,3:end-2),medR(4,3:end),'-y','LineWidth',2)
plot(med(n,3:end-2),medR(5,3:end),'-y','LineWidth',2)
plot(med(p,3:end-2),medR(6,3:end),'-y','LineWidth',2)
plot(cor(f,3:end-2),(7*ones(size(corR(5,3:end)))),'--k','LineWidth',2)
hold off

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

[q,r]=size(cor);
cor1 = -22*ones(q,1);
cor = horzcat(cor1,cor);
[q,r]=size(cor);
timestep = 12.5;%seconds
tminut = timestep/60;%minute



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

for j=1:50
     corr=(cor(b,j+1)-cor(b,j))/12.5;
     corR(1,j) = corr;

     corr=(cor(c,j+1)-cor(c,j))/12.5;
     corR(2,j) = corr;
     
     corr=(cor(d,j+1)-cor(d,j))/12.5;
     corR(3,j) = corr;
     
     corr=(cor(e,j+1)-cor(e,j))/12.5;
     corR(4,j) = corr;

     corr=(cor(f,j+1)-cor(f,j))/12.5;
     corR(5,j) = corr;
     
     corr=(cor(k,j+1)-cor(k,j))/12.5;
     corR(6,j) = corr;
end
for j=51:r-2
    
    psb=horzcat(cor(b,j-2),cor(b,j-1),cor(b,j),cor(b,j+1),cor(b,j+2));
    tsb=time(j-2:j+2);
    sb = polyfit(tsb,psb,1);
    corR(1,j) = sb(1);
    
    psc=horzcat(cor(c,j-2),cor(c,j-1),cor(c,j),cor(c,j+1),cor(c,j+2));
    sc = polyfit(tsb,psc,1);
    corR(2,j) = sc(1);
    
    psd=horzcat(cor(d,j-2),cor(d,j-1),cor(d,j),cor(d,j+1),cor(d,j+2));
    sd = polyfit(tsb,psd,1);
    corR(3,j) = sd(1);
    
    pse=horzcat(cor(e,j-2),cor(e,j-1),cor(e,j),cor(e,j+1),cor(e,j+2));
    se = polyfit(tsb,pse,1);
    corR(4,j) = se(1);
    
    psf=horzcat(cor(f,j-2),cor(f,j-1),cor(f,j),cor(f,j+1),cor(f,j+2));
    sf = polyfit(tsb,psf,1);
    corR(5,j) = sf(1);
    
    psk= horzcat(cor(k,j-2),cor(k,j-1),cor(k,j),cor(k,j+1),cor(k,j+2));
    sk = polyfit(tsb,psk,1);
    corR(6,j) = sk(1);
    
end

[q,r]=size(med);
med1 = -22*ones(q,1);
med = horzcat(med1,med);
[q,r]=size(med);

for j=3:r-2
    psa=horzcat(med(a,j-2),med(a,j-1),med(a,j),med(a,j+1),med(a,j+2));
    sa = polyfit(tsb,psa,1);
    medR(1,j) = sa(1);
    
    psg=horzcat(med(g,j-2),med(g,j-1),med(g,j),med(g,j+1),med(g,j+2));
    sg = polyfit(tsb,psg,1);
    medR(2,j) = sg(1);
    
    psh=horzcat(med(h,j-2),med(h,j-1),med(h,j),med(h,j+1),med(h,j+2));
    sh = polyfit(tsb,psh,1);
    medR(3,j) = sh(1);
    
    psm=horzcat(med(m,j-2),med(m,j-1),med(m,j),med(m,j+1),med(m,j+2));
    sm = polyfit(tsb,psm,1);
    medR(4,j) = sm(1);
    
    psn=horzcat(med(n,j-2),med(n,j-1),med(n,j),med(n,j+1),med(n,j+2));
    sn = polyfit(tsb,psn,1);
    medR(5,j) = sn(1);
    
    psp=horzcat(med(p,j-2),med(p,j-1),med(p,j),med(p,j+1),med(p,j+2));
    sp = polyfit(tsb,psp,1);
    medR(6,j) = sp(1);
end
%%


%chunks instead of individual lines
%
corR=-corR*60;
medR=-medR*60;
figure(3)
hold on
plot(cor(b,1:end-2),corR(1,1:end),'-c','LineWidth',2)
hold on
plot(cor(c,1:end-2),corR(2,1:end),'-c','LineWidth',2)
plot(cor(d,1:end-2),corR(3,1:end),'-c','LineWidth',2)
plot(cor(e,1:end-2),corR(4,1:end),'-c','LineWidth',2)
plot(cor(f,1:end-2),corR(5,1:end),'-c','LineWidth',2)
plot(cor(k,1:end-2),corR(6,1:end),'-c','LineWidth',2)
%plot(cor(f,1:end-2),(0.1*ones(size(corR(5,3:end)))),'--k','LineWidth',2)



figure(3)
hold on
plot(med(a,3:end-2),medR(1,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(time(3:end-2),medR(2,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(g,3:end-2),medR(2,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(h,3:end-2),medR(3,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(m,3:end-2),medR(4,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(n,3:end-2),medR(5,3:end),'-','color',[1 0.5 0],'LineWidth',2)
plot(med(p,3:end-2),medR(6,3:end),'-','color',[1 0.5 0],'LineWidth',2)
%plot(cor(f,3:end-2),(7*ones(size(corR(5,3:end)))),'--k','LineWidth',2)
hold off

