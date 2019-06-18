load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_8_21_t1000_cortex.mat

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\cortexx.csv'
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\cortexy.csv'
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\cortexz.csv'
Z = importdata(filename);

Ccoords = horzcat(X,Y,Z);

[m,n]=size(corM);

%scatter3(X,Y,Z)

timestep = 100;%seconds
tminut = timestep/60;%minutes

for j=1:n-1
    curt = corM(:,j);
    next = corM(:,j+1);
    corR(:,j) = (curt-next)/tminut;
end

[k,l]=size(corR);
colorcor = zeros(k,l);
for i=1:k
    for h =1:l
        if corR(i,h)<=0.1,
            colorcor (i,h) = 0;
        else
            colorcor(i,h) = 1;
        end
    end
end
map = [1,0,0 
       0,0,1];

figure
subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);
scatter3(X,Y,Z,2,colorcor(:,2))%,'filled','markeredgecolor','k')
%view(40,75)
colormap(map)

loops = 109
time = 100*linspace(1,109,109);
F(loops) = struct('cdata',[],'colormap',[]);
figure
subplot(1,2,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);

for i=1:loops
    tmax(i) = max(corM(:,i))
    tmin(i) = min(corM(:,i))
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
figure
for j =1:loops
     subplot(1,2,1)
     set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 15, 9]);
    scatter3(X,Y,Z,2,corM(:,j))
    colormap(jet)
    colorbar
    subplot(1,2,2)
    plot(time,tmax,'-r')
    hold on
    plot(time,tmin,'-b')
    plot(time(j),tmax(j),'o','markerfacecolor','r','markeredgecolor','k')
    plot(time(j),tmin(j),'o','markerfacecolor','b','markeredgecolor','k')
    hold off
    drawnow
    G(j) = getframe(gcf);
end
fig = figure;
%subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 15, 9]);
movie(fig,F,2)

fig = figure;
%subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 15, 9]);
movie(fig,G,2)

time = 100*linspace(1,109,109);
loops = 109
F(loops) = struct('cdata',[],'colormap',[]);







