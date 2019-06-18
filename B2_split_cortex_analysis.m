clear all
clc

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
plot(ccoord_slice(:,2),ccoord_slice(:,3),'.r')
hold on
plot(mcoord_slice(:,2),mcoord_slice(:,3),'.b')

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
    

[k,l]=size(corR);
colorcor = zeros(k,l);
for i=1:k
    for h =1:l
        if (corR(i,h)<=0.1) && (cor(i,h)<=-40) && (cor(i,h)>=-80)
            colorcor(i,h) = 0;
        else
            colorcor(i,h) = 1;
        end
    end
end
colorc=colorcor;
mapc = [1,0,0 
       0,0,1];

[k,l]=size(medR);
colormed = zeros(k,l);
for i=1:k
    for h =1:l
        if (medR(i,h)<=7) && (med(i,h)<=-40) && (med(i,h)>=-80)
            colormed(i,h) = 0;
        else
            colormed(i,h) = 1;
        end
    end
end
colorm=colormed;
mapm = [1,0,0 
       0,0,1];


    
    
    
    
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

% figure
% a=8958;
% b=79;
% c=24429;
% d= 10251;
% e=46159;
% f=781;
% g=309;
% h=10839;
% i=25000;
% j=272;
figure

scatter3(Ccoords(:,1),Ccoords(:,2),Ccoords(:,3),2,'.k')
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);
hold on
scatter3(Ccoords(b,1),Ccoords(b,2),Ccoords(b,3),100,'h','markerfacecolor','c','markeredgecolor','c')
scatter3(Ccoords(c,1),Ccoords(c,2),Ccoords(c,3),100,'h','markerfacecolor','b','markeredgecolor','b')
scatter3(Ccoords(d,1),Ccoords(d,2),Ccoords(d,3),100,'h','markerfacecolor',[.49 1 .63],'markeredgecolor',[.49 1 .63])
scatter3(Ccoords(e,1),Ccoords(e,2),Ccoords(e,3),100,'h','markerfacecolor',[0.5 0.5 1],'markeredgecolor',[0.5 0.5 1])
scatter3(Ccoords(f,1),Ccoords(f,2),Ccoords(f,3),100,'h','markerfacecolor', 'g','markeredgecolor','g')
scatter3(Ccoords(k,1),Ccoords(k,2),Ccoords(k,3),100,'h','markerfacecolor',[0.5 0 1],'markeredgecolor',[0.5 0 1])
scatter3(Mcoords(a,1),Mcoords(a,2),Mcoords(a,3),100,'h','markerfacecolor','r','markeredgecolor','r')
scatter3(Mcoords(g,1),Mcoords(g,2),Mcoords(g,3),100,'h','markerfacecolor',[0.5 0 0],'markeredgecolor',[0.5 0 0])
scatter3(Mcoords(h,1),Mcoords(h,2),Mcoords(h,3),100,'h','markerfacecolor',[1 0.5 0],'markeredgecolor',[1 0.5 0])
scatter3(Mcoords(m,1),Mcoords(m,2),Mcoords(m,3),100,'h','markerfacecolor',[1 0.5 0.5],'markeredgecolor',[1 0.5 0.5])
scatter3(Mcoords(n,1),Mcoords(n,2),Mcoords(n,3),100,'h','markerfacecolor','m','markeredgecolor','m')
scatter3(Mcoords(p,1),Mcoords(p,2),Mcoords(p,3),100,'h','markerfacecolor','y','markeredgecolor','y')

%light blue green [.49 1 .63]
%blue purple [0.5 0.5 1]
%purple [0.5 0 1]

figure
plot(time,cor(b,:),'-c','LineWidth',2)
hold on
plot(time,cor(c,:),'-b','LineWidth',2)
plot(time,cor(d,:),'MarkerEdgeColor',[0.49 1 0.63],'LineWidth',2)
plot(time,cor(e,:),'markeredgecolor',[0.5 0.5 1],'LineWidth',2)
plot(time,cor(f,:),'-g','LineWidth',2)
plot(time,cor(k,:),'MarkerEdgeColor',[0.5 0 1],'LineWidth',2)
plot(time,med(a,:),'-r','LineWidth',2)
plot(time,med(g,:),'MarkerEdgeColor',[0.5 0 0],'LineWidth',2)
plot(time,med(h,:),'MarkerEdgeColor',[1 0.5 0],'LineWidth',2)
plot(time,med(m,:),'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',2)
plot(time,med(n,:),'-m','LineWidth',2)
plot(time,med(p,:),'-y','LineWidth',2)

figure
plot(cor(b,1:end-1),corR(b,:),'-c','LineWidth',2)
hold on
plot(cor(c,1:end-1),corR(c,:),'-b','LineWidth',2)
plot(cor(d,1:end-1),corR(d,:),'MarkerEdgeColor',[0.49 1 0.63],'LineWidth',2)
plot(cor(e,1:end-1),corR(e,:),'markeredgecolor',[0.5 0.5 1],'LineWidth',2)
plot(cor(f,1:end-1),corR(f,:),'-g','LineWidth',2)
plot(cor(k,1:end-1),corR(k,:),'MarkerEdgeColor',[0.5 0 1],'LineWidth',2)
plot(med(a,1:end-1),medR(a,:),'-r','LineWidth',2)
plot(med(g,1:end-1),medR(g,:),'MarkerEdgeColor',[0.5 0 0],'LineWidth',2)
plot(med(h,1:end-1),medR(h,:),'MarkerEdgeColor',[1 0.5 0],'LineWidth',2)
plot(med(m,1:end-1),medR(m,:),'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',2)
plot(med(n,1:end-1),medR(n,:),'-m','LineWidth',2)
plot(med(p,1:end-1),medR(p,:),'-y','LineWidth',2)


figure
plot(cor(b,2:end),corR(b,:),'-c','LineWidth',2)
hold on
plot(cor(c,2:end),corR(c,:),'-b','LineWidth',2)
plot(cor(d,2:end),corR(d,:),'MarkerEdgeColor',[0.49 1 0.63],'LineWidth',2)
plot(cor(e,2:end),corR(e,:),'markeredgecolor',[0.5 0.5 1],'LineWidth',2)
plot(cor(f,2:end),corR(f,:),'-g','LineWidth',2)
plot(cor(k,2:end),corR(k,:),'MarkerEdgeColor',[0.5 0 1],'LineWidth',2)
plot(cor(f,2:end),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)

figure
hold on
plot(med(a,2:end),medR(a,:),'-r','LineWidth',2)
plot(med(g,2:end),medR(g,:),'MarkerEdgeColor',[0.49 0 0],'LineWidth',2)
plot(med(h,2:end),medR(h,:),'MarkerEdgeColor',[1 0.5 0],'LineWidth',2)
plot(med(m,2:end),medR(m,:),'MarkerEdgeColor',[1 0.5 0.5],'LineWidth',2)
plot(med(n,2:end),medR(n,:),'-m','LineWidth',2)
plot(med(p,2:end),medR(p,:),'-y','LineWidth',2)
plot(cor(f,2:end),(7*ones(size(corR(f,:)))),'--k','LineWidth',2)


% 
% figure
% plot(cor(a,:),corR(a,:),'-m','LineWidth',2)
% hold on
% plot(cor(h,:),corR(h,:),'-b','LineWidth',2)
% plot(cor(i,:),corR(i,:),'-c','LineWidth',2)
% plot(cor(d,:),corR(d,:),'-g','LineWidth',2)
% plot(cor(e,:),corR(e,:),'-r','LineWidth',2)
% plot(cor(f,:),corR(f,:),'-k','LineWidth',2)
% plot(cor(g,:),corR(g,:),'MarkerEdgeColor',[.49 1 .63],'LineWidth',2)
% plot(cor(j,:),corR(j,:),'MarkerEdgeColor',[.5 0.5 0],'LineWidth',2)
% plot(cor(f,:),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)


