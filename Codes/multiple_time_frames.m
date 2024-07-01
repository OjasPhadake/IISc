% Code Summary:
% Initial code in which I tried to read each line one by one and then
% process it. Later learnt that importdata function does the same work much
% more efficiently. Spent a lot of time on this, as I grappled with
% fighting the time format from string to char to int. 

%% Clearing previous data
clc;
clear all;
close all;

%% 


vz_arr = zeros(12, 1);
j=1;

for time = 63500:500:69000
    A = fileread(convertStringsToChars("C:\IISc\DEM\1case\post\particles_" + time + ".liggghts"));
    data = splitlines(A);
    data = data(10:end);

    for i = 1:length(data)-1
        line = convertCharsToStrings(data{i});
        all_lines(i, :) = str2num(line);
    end
    
    vz = all_lines(:, 8);
    vz_arr(j, 1) = sum(vz)/length(data);
    j = j + 1;
end

plot(vz_arr)
