% Code Summary:
% Short code for reading time series data. Can use this as modular to read
% more data with time. Simply plotting granular temperature wrt time.


%% Clearing previous data
clc;
close all;
clear all;

%% Load the data

filenums = 100000:1000:150000;
gt = zeros(size(filenums, 2), 1);
i = 1;

for time=filenums
    filename = "C:\IISc\DEM\1case\post\particles_" + time + ".liggghts";
    file = importdata(filename, " ", 9);
    data = file.data;
    
    vx = data(:, 8);
    vy = data(:, 9);
    vz = data(:, 10);

    vxm = vx - mean(vx);
    vym = vy - mean(vy);
    vzm = vz - mean(vz);

    Tgran = 1/3 * mean((vxm.^2 + vym.^2 +vzm.^2));
    gt(i) = Tgran;
    i = i + 1;
end

%% Plotting 

figure
plot(filenums, gt)
xlabel("File number")
ylabel("Granular Temperature")
title("Granular temp vs file number")
