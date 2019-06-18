load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_8_15_t8000_medulla

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\medx.csv'
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\medy.csv'
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\medz.csv'
Z = importdata(filename);

Ccoords = horzcat(X,Y,Z);

[m,n]=size(medM);

%scatter3(X,Y,Z)

timestep = 100;%seconds
tminut = timestep/60;%minutes

for j=1:n-1
    curt = medM(:,j);
    next = medM(:,j+1);
    medR(:,j) = (curt-next)/tminut;
end

[k,l]=size(medR);
colormed = zeros(k,l);
for i=1:k
    for h =1:l
        if medR(i,h)<=2,
            colormed (i,h) = 0;
        else
            colormed(i,h) = 1;
        end
    end
end
map = [1,0,0 
       0,0,1];

figure
subplot(1,3,1)
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 6.83, 9]);
scatter3(X(76790:end),Y(76790:end),Z(76790:end),2,colormed(:,2))%,'filled','markeredgecolor','k')
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
    tmax(i) = max(medM(:,i))
    tmin(i) = min(medM(:,i))
end
for j = 1:loops
    subplot(1,2,1)
    set(gcf, 'Units', 'inches');
    set(gcf, 'Position', [0, 0, 15, 9]);
    colormap(map)
    scatter3(X(76790:end),Y(76790:end),Z(76790:end),3,colormed(:,j))%,'filled','markeredgecolor','k')
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