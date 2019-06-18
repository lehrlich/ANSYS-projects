%load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_8_21_t1000_cortex.mat

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medonex.csv';
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medoney.csv';
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\medonez.csv';
Z = importdata(filename);

Ccoords = horzcat(X,Y,Z);
%Ccoords = Ccoords(63206:end,:);


CPA_contents_centr = [1.0707e-5 -7.2976e-6 6.3525e-2];
kidney_centroid = [-5.5334e-3 3.81e-3 4.9603e-2];
CPA = repmat(CPA_contents_centr,length(Ccoords),1);
cent = repmat(kidney_centroid,length(Ccoords),1);

dist = sqrt(((Ccoords(:,1)-CPA(:,1)).^2)+((Ccoords(:,2)-CPA(:,2)).^2)+((Ccoords(:,3)-CPA(:,3)).^2));
dist_c = sqrt(((Ccoords(:,1)-cent(:,1)).^2)+((Ccoords(:,2)-cent(:,2)).^2)+((Ccoords(:,3)-cent(:,3)).^2));
sliceM = zeros(size(Ccoords));
for i=1:length(Ccoords)
    if Ccoords(i,1)>=-5.334e-3 & Ccoords(i,1)<=-5.0694e-3
        sliceM(i,:) = 1;
%     elseif Ccoords(i,1)>= -5.2899141e-3;
%         sliceM(i,:) = 1;
    else
        sliceM(i,:) = NaN;
    end
end

dist_slice = dist.*sliceM(:,1);
dist_c_slice = dist_c.*sliceM(:,1);
coord_slice = Ccoords.*sliceM;

[min_dist,I_7] = min(dist_slice)
[max_dist,I_10]=max(dist_slice)
[min_dist,I_1]=min(dist_c_slice)
plot(coord_slice(:,2),coord_slice(:,3),'.r')
