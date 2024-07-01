% Code Summary:
% I was bored, and hence I decided to write a toy code to calculate the
% repose angle of a given heap. So, I've done radial and z binning and
% whereas the z binning was accurate, the r binning wasn't giving me the
% desired results, so I've tweaked the work for a (currently) seeming
% arbitrary criterion for num_particles(i) and num_particles(i+1) so that
% it'll detect if there is a change. One nice thing which I did however,
% was to remove the particles which are in the lowermost layer. I did this
% by find(z(A) < 1.2*radius). This allowed me to remove the single
% molecules which are lying here and there, and made the heap much smaller,
% and better. 

%% Clearing previous data
clc;
close all;
clear variables;

%% Load the file

rep_arr = [];
for time=361000:1000:390000
str = "post\particles_" + time + ".liggghts";
file = importdata(str, " ", 9);
data = file.data;
clear file;

X = data(:, 3);
Y = data(:, 4);
Z = data(:, 5);
vz = data(:, 8);
radius = data(1, end-1);
dp = 2*radius;

%% Radial stuff
zbool = (Z<0.036);
x = X(zbool); y = Y(zbool); z = Z(zbool);
r = sqrt(x.^2 + y.^2);

radialmesh = linspace(0, max(z), 100); % Dividing the radius into 100 small parts
zmesh = linspace(0, max(z), 100); % Dividing z
num_particles = zeros(1, 100); % Radially

for i=1:length(num_particles) - 1
    A = find(r > radialmesh(i) & r < radialmesh(i+1));
    num_particles(i) = length(A);
    B = find(z(A) < 1.2*radius);
    num_particles(i) = num_particles(i) - length(B);
end

for i=1:length(num_particles)
    if ((num_particles(i) - num_particles(i+1))/num_particles(i+1) > 5 ) 
        width = radialmesh(i);
        break
    end
end
% For z heights, I'll try looking at z velocity. 

num_particles = zeros(1, 100); 

zvels = zeros(1, 100);
for i=1:length(num_particles) - 1
    rows = find(z > zmesh(i) & z < zmesh(i+ 1));
    num_particles(i) = length(rows);
    zvels(i) = sum(vz(rows));
end

normzvels = zvels./num_particles;
height = zmesh(normzvels > min(z(z>=0)));
if (isempty(height))
  continue;
else
  height = height(1);
end

repose = atan(height/width)*180/pi;
rep_arr = [rep_arr, repose];
end

median(rep_arr)
mode(rep_arr)
mean(rep_arr)
