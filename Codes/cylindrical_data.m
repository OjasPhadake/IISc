%% Task 2 of 2D cylindrical flow - Converting the C++ code into MATLAB code
clc;
close all;
clear all;

%% Load the data

% ntimes = 10;
% data = zeros(685, 20, ntimes);
% hlines = 9;
% fid = fopen("dump_m0.18_v-1.61.atom");
% 
% for i=1:ntimes
%     if feof(fid)
%         break;
%     end
%     file = importdata("dump_m0.18_v-1.61.atom", " ", hlines);
%     hlines = hlines + 694;
%     data(:, :, i) = file.data;
% end

% The above code is for reading the dump file with large size. It is
% basically a fix file which is getting outputted in the same file instead
% of multiple files. 

% However, as the number of atoms will keep on changing, hence what I
% should do it, importdata(filename, " ", 3) which will give me the number
% of atoms, and then for the next timestep, add the numberofatoms + 9 of
% the previous step to start reading the data for the next time step
% PATHETIC ATTEMPT at copying a C++ code to MATLAB, which went
% unsuccessful. It'd have taken much effort to actually transfer it in
% exact copy. However, the generic_code I think does justice to the spirit
% of inquiry and creates bins

%% Defining variables

id=0; 
x=0; y=0; z=0; 
ux=0; uy=0; uz=0; us=0;
vx=0; vy=0; vz=0; 
m=0; s=0; s_old=0; 
fx=0; fy=0; fz=0; 
omx=0; omy=0; omz=0; 
cxx=0; cyy=0; czz=0; 
cxy=0; cyz=0; cxz=0;
r=0; count=0;

xb=0.001; yb=0.0101; zb=0.001; % dimensions of the bin
xl=-0.08; yl=-0.0001; zl=-0.15;
xu=0.08; yu=0.01; zu=0.10;

tstart = 100000; tstop = 300000; Vc=0; rho=2500;
tbx = 0; tby = 0; tbz = 0; tbin = 0;

%% Code begins

tbx = gif((xu - xl), xb);
tby = gif((yu - yl), yb);
tbz = gif((zu - zl), zb);

tbin = tbx*tby*tbz;

fprintf("The values for bin params are %0.1d, %0.1d, %0.1d, %0.1d\n", tbin, tbx, tby, tbz)
Vc = xb*yb*zb;

% Read from the file dump_m0.115_v-1.61.atom
% Write to the file gstressout_m0.16_v2.2.dat

Vbin = zeros(tbin, 7, 200);
inf2 = zeros(tbin, 2);

timestep = 0; No_of_atoms = 0; flag = 1;
hlines = 9;
fid = fopen("dump_m0.18_v-1.61.atom");
data = zeros(685, 20);

while (flag)
    Vbin = zeros(tbin, 7, 200);

    t = importdata("dump_m0.18_v-1.61.atom");
    timestep = t.data;

    no = importdata("dump_m0.18_v-1.61.atom", " ", 3);
    No_of_atoms = no.data;

    if feof(fid)
        break;
    end

    file = importdata("dump_m0.18_v-1.61.atom", " ", hlines);
    hlines = hlines + No_of_atoms + 9;
    data(:, :) = file.data;
    
    % Now, I have the data matrix named data which contains : as columns
    % ATOMS id type x y z v_tbin vx vy vz v_s fx fy fz c_1[1] c_1[2] c_1[3] c_1[4] c_1[5] c_1[6] radius 
    
    if (timestep>tstop)
        flag = 0;
    elseif (timestep<tstart)
        % How to implement: Ignore current line?
        % inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        continue; % Is this the right way to proceed? 
    else
        count = count + 1; % no of timesteps processed
        x = data(:, 3);
        y = data(:, 4);
        z = data(:, 5);
        vx = data(:, 7);
        vy = data(:, 8);
        vz = data(:, 9);

        cxx = data(:, 14);
        cyy = data(:, 15);
        czz = data(:, 16);

        cxy = data(:, 17);
        cxz = data(:, 18);
        cyz = data(:, 19);
        r = data(:, 20);

        % bin = ceil((x - xl) / xb) + floor((y - yl) / yb) * tbx + floor((z - zl) / zb) * tbx * tby; 
        
    end


end


%% Helper functions

function [op] = gif(x, y)
a = x/y;
b = floor(a);

if b*y ~= x
    op = b + 1;
else
    op = floor(x/y);
end
end
