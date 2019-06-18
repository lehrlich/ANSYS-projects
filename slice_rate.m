clear all
clc

load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_8_21_12p5_cortex.mat

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexx.csv'
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexy.csv'
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexz.csv'
Z = importdata(filename);


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

corR = horzcat(corR{1},corR{2},corR{3},corR{4},corR{5},corR{6},corR{7},corR{8});
%plot cooling rates along a slice of kidney
%along x = 0
j=1;
V=zeros(size(corR));
C=zeros(size(X,3));
for i = 1:length(X)
    if abs(X(i)) <=   9.0970e-03

        V(j,:) = corR(i,:);
        C(j,1) = X(i);
        C(j,2) = Y(i);
        C(j,3)=Z(i);
        j=j+1;
    end
end

plot(C(:,1),C(:,2),'o')
