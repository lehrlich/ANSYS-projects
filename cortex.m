clear all
clc

for g = 1:4
    
    filename = strcat('C:\Users\Student\Desktop\3omega\ANSYS_kidney\parse_results\A2_split_12p5_results\med2result',num2str(g),'.csv');
    %('C:\Users\Student\Desktop\kidney_ANSYS\kidney_B2_8_24\kidney_B2_10_11_t8000_12p5_split\kidney_B2_10_11_t8000_12p5_split_files\user_files\med1result',num2str(g),'.csv');
    %strcat('C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_10_11_t11000_12p5_split\kidney_A2_10_11_t11000_12p5_split_files\user_files\cortex1result',num2str(g),'.csv');
    %filename = ('C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_10_11_t11000_12p5_split\kidney_A2_10_11_t11000_12p5_split_files\user_files\med1result1.csv')
    %'C:\Users\Student\Desktop\kidney_ANSYS\kidney_rA2_8_21\kidney_rA2_8_21_t1000_files\user_files\cortexresult.csv';
    %'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t11000_t12p5s\kidney_A2_8_15_t11000_t12p5s_files\user_files\cortexresult.csv';
    %'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\cortexresult.csv'
    M = importdata(filename);
    T = struct2cell(M);
    T2 = T{2,1};
    %T2c = cell2struct(T2);
    i=1;
    j=1;
    n=1;

    j=0;
    n=1;
    [x,z]= size(T2);
    medM =zeros(3065,z);
for n = 1:length(T2);
    medT = cell2mat(T2(n));
    Q=str2num(medT);
    if isempty(Q)==1;
        j=j+1;
        i=1;
        %wiret(j)=str2num(cell2mat(T2(n+1)));
    elseif Q ==0
    else medM(i,j) = Q;
        i=i+1;
    end
end

medM( ~any(medM,2), : ) = [];  
med{g}=medfahyM;
end

filename = 'A2_split_12p5_med2';
save(filename,'med');