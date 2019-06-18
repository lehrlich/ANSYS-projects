load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_8_21_12p5_cortex.mat

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexx.csv'
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexy.csv'
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexz.csv'
Z = importdata(filename);

%X =X(1:52964,1);

Ccoords = horzcat(X,Y,Z);

for c = 1:8;
[m,n]=size(cor{c});

%scatter3(X,Y,Z)

timestep = 12.5;%seconds
tminut = timestep/60;%minutes

for j=1:n-1
    curt = cor{c}(:,j);
    next = cor{c}(:,j+1);
    corR{c}(:,j) = abs((curt-next)/tminut);
end

[k,l]=size(corR{c});
colorcor = zeros(k,l);
for i=1:k
    for h =1:l
        if (corR{c}(i,h)<=0.1) && (cor{c}(i,h)<=-40) && (cor{c}(i,h)>=-80)
            colorcor(i,h) = 0;
        else
            colorcor(i,h) = 1;
        end
    end
end
colorc{c}=colorcor;
map = [1,0,0 
       0,0,1];
end

colorcor = horzcat(colorc{1},colorc{2},colorc{3},colorc{4},colorc{5},colorc{6},colorc{7},colorc{8});
cor = horzcat(cor{1},cor{2},cor{3},cor{4},cor{5},cor{6},cor{7},cor{8});
corR = horzcat(corR{1},corR{2},corR{3},corR{4},corR{5},corR{6},corR{7},corR{8});
figure
subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);
scatter3(X,Y,Z,2,colorcor(:,2))%,'filled','markeredgecolor','k')
%view(40,75)
colormap(map);

loops = 872;
time = 12.5*linspace(1,872,872);
F(loops) = struct('cdata',[],'colormap',[]);
figure
subplot(1,2,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);

for i=1:loops
    tmax(i) = max(cor(:,i));
    tmin(i) = min(cor(:,i));
end
for j = 1:loops
    subplot(1,2,1)
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 15, 9]);
    colormap(map)
    scatter3(X,Y,Z,2,colorcor(:,j))%,'filled','markeredgecolor','k')
    subplot(1,2,2)
    plot(time,tmax,'-r')
    hold on
    plot(time,tmin,'-b')
    plot(time(j),tmax(j),'o','markerfacecolor','r','markeredgecolor','k')
    plot(time(j),tmin(j),'o','markerfacecolor','b','markeredgecolor','k')
    hold off
    drawnow
    F(j) = getframe(gcf);
end
%figure
% for j =1:loops
%      subplot(1,2,1)
%      set(gcf, 'Units', 'inches');
%     set(gcf, 'Position', [0, 0, 15, 9]);
%     scatter3(X,Y,Z,2,corM(:,j))
%     colormap(jet)
%     colorbar
%     subplot(1,2,2)
%     plot(time,tmax,'-r')
%     hold on
%     plot(time,tmin,'-b')
%     plot(time(j),tmax(j),'o','markerfacecolor','r','markeredgecolor','k')
%     plot(time(j),tmin(j),'o','markerfacecolor','b','markeredgecolor','k')
%     hold off
%     drawnow
%     G(j) = getframe(gcf);
% end
fig = figure;
%subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 15, 9]);
movie(fig,F,2)
% 
% fig = figure;
% %subplot(1,3,1)
% set(gcf, 'Units', 'inches');
% set(gcf, 'Position', [0, 0, 15, 9]);
% movie(fig,G,2)
% % 
% time = 100*linspace(1,109,109);
% loops = 109
% F(loops) = struct('cdata',[],'colormap',[]);
% 
time = 12.5*linspace(1,880,880);


figure
a=8958;
b=79;
c=24429;
d= 10251;
e=46159;
f=781;
g=309;
h=10839;
i=25000;
j=272;

scatter3(X,Y,Z,2,'.k')
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);
hold on
scatter3(X(a),Y(a),Z(a),100,'h','markerfacecolor','m','markeredgecolor','m')
scatter3(X(h),Y(h),Z(h),100,'h','markerfacecolor','b','markeredgecolor','b')
scatter3(X(i),Y(i),Z(i),100,'h','markerfacecolor','c','markeredgecolor','c')
scatter3(X(d),Y(d),Z(d),100,'h','markerfacecolor','g','markeredgecolor','g')
scatter3(X(e),Y(e),Z(e),100,'h','markerfacecolor','r','markeredgecolor','r')
scatter3(X(f),Y(f),Z(f),100,'h','markerfacecolor','k','markeredgecolor','k')
scatter3(X(g),Y(g),Z(g),100,'h','markerfacecolor',[.49 1 .63],'markeredgecolor',[.49 1 .63])
scatter3(X(j),Y(j),Z(j),100,'h','markerfacecolor',[.5 0.5 0],'markeredgecolor',[.5 0.5 0])

figure
plot(time,cor(a,:),'-m','LineWidth',2)
hold on
plot(time,cor(h,:),'-b','LineWidth',2)
plot(time,cor(i,:),'-c','LineWidth',2)
plot(time,cor(d,:),'-g','LineWidth',2)
plot(time,cor(e,:),'-r','LineWidth',2)
plot(time,cor(f,:),'-k','LineWidth',2)
plot(time,cor(g,:),'MarkerEdgeColor',[.49 1 .63],'LineWidth',2)
plot(time,cor(j,:),'MarkerEdgeColor',[.5 0.5 0],'LineWidth',2)

cor=horzcat(cor(:,1:109),cor(:,110:219),cor(:,220:328),cor(:,329:437),cor(:,438:546),cor(:,547:655),cor(:,656:764),cor(:,765:872));

figure
plot(cor(a,:),corR(a,:),'-m','LineWidth',2)
hold on
plot(cor(h,:),corR(h,:),'-b','LineWidth',2)
plot(cor(i,:),corR(i,:),'-c','LineWidth',2)
plot(cor(d,:),corR(d,:),'-g','LineWidth',2)
plot(cor(e,:),corR(e,:),'-r','LineWidth',2)
plot(cor(f,:),corR(f,:),'-k','LineWidth',2)
plot(cor(g,:),corR(g,:),'MarkerEdgeColor',[.49 1 .63],'LineWidth',2)
plot(cor(j,:),corR(j,:),'MarkerEdgeColor',[.5 0.5 0],'LineWidth',2)
plot(cor(f,:),(0.1*ones(size(corR(f,:)))),'--k','LineWidth',2)
