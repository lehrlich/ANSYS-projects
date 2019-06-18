clear all
clc

divis = 77;
filename = 'C:\Users\Student\Desktop\kidney_ANSYS\kidney_A2_8_15\kidney_A2_8_15_t8000_files\user_files\medresult.csv'
%'E:\3Omega\ANSYS_cuvette_analysis\simulations\v7x1_square\v7x1_square_files\user_files\pathresult.csv'
%'E:\3Omega\ANSYS_cuvette_analysis\simulations\Thermal_model_CPA_cuvette_v4_const_timestep\Thermal_model_CPA_cuvette_v4_const_timestep_files\user_files\pathresult.csv';
%'E:\3Omega\ANSYS_cuvette_analysis\simulations\Thermal_model_CPA_cuvette_v4_without_hgen\Thermal_model_CPA_cuvette_v4_without_hgen_files\user_files\pathresult.csv';
%'E:\3Omega\ANSYS_cuvette_analysis\simulations\Thermal_model_CPA_cuvette_v4\Thermal_model__files\user_files\pathresult.csv';
% 'E:\3Omega\ANSYS_cuvette_analysis\Thermal_model_CPA_cuvette_v3_files\user_files\pathresult.csv';
M = importdata(filename);
T = struct2cell(M);
T2 = T{2,1};
%T2c = cell2struct(T2);
i=1;
j=1;
n=1;
%%
%%

%%
%%
j=0;
n=1;

for n = 1:length(T2);
    medT = cell2mat(T2(n));
    Q=str2num(medT);
    if isempty(Q)==1;
        j=j+1;
        i=1;
        %wiret(j)=str2num(cell2mat(T2(n+1)));
    else medM(i,j) = Q;
        i=i+1;
    end
end

medM( ~any(medM,2), : ) = [];  %rows
med{g}=medM
