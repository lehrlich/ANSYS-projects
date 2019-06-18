%load C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_8_21_t1000_cortex.mat

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonex.csv';
X = importdata(filename);

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideoney.csv';
Y = importdata(filename);

filename = 'C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\rA2_split_12p5_results\corsideonez.csv';
Z = importdata(filename);

Ccoords = horzcat(X,Y,Z);
%Ccoords = Ccoords(63206:end,:);


CPA_contents_centr = [-8.8619e-8 -1.3684e-7 0.022859];
kidney_centroid = [-0.0019897 0.001372 0.017856];
CPA = repmat(CPA_contents_centr,length(Ccoords),1);
cent = repmat(kidney_centroid,length(Ccoords),1);

dist = sqrt(((Ccoords(:,1)-CPA(:,1)).^2)+((Ccoords(:,2)-CPA(:,2)).^2)+((Ccoords(:,3)-CPA(:,3)).^2));
dist_c = sqrt(((Ccoords(:,1)-cent(:,1)).^2)+((Ccoords(:,2)-cent(:,2)).^2)+((Ccoords(:,3)-cent(:,3)).^2));
sliceM = zeros(size(Ccoords));
for i=1:length(Ccoords)
    if Ccoords(i,1)>=-1.9203e-3 & Ccoords(i,1)<=-1.7e-3
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
%[max_dist,I_10]=max(dist_slice)
[min_dist,I_1]=min(dist_c_slice)
%[max_dist,I_2]=max(Ccoords(:,3))
%[max_dist,I_3]=max(Ccoords(:,2))
%[min_dist,I_5]=min(Z(Z>0))
plot(-1*coord_slice(:,2),coord_slice(:,3),'.b')
