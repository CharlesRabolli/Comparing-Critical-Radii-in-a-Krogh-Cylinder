%% Initial

clear all;
clc;

D = 1*10^(-7); %cm^2/s Diffusion Coefficient
V = 0.002; %cm/s Velocity
rc = 0.0005; %cm capillary radius
tm = 5*10^(-5); %cm capillary wall thickness
Ko = 5.75*10^(-5); %cm/s Overall Mass Transfer Rate
Ro = 0.01; %umole/cm^3 * s Metabolic Consumption Rate
Co = 5; %umole/cm^3
t = [0.001:0.01:0.101]; %varying z from 0.001 to 0.1 by 0.01

r_crit = @(y) y.*(rc + tm);
Krogh = @(y,z) y.^2.*log(y.^2) - y.^2 + 1 - ((4*D*Co)/(Ro*(rc+tm)^2)) + 4*(D/(V*rc^2)).*(y.^2-1)*z+(2*D/(rc*Ko)).*(y.^2-1);
dKdy = @(y,z) 2.*y+2.*y.*log(y.^2) - 2*y +8*(D/(V*rc^2)).*y.*z + (4*D/(rc*Ko)).*y;

tol = 10^(-7);
maxit=100;
for i=1:11
    z=t(i);
    y=0.001; q=1;
    while q<100
        ynew = y - (Krogh(y,z)/dKdy(y,z));
        Error = abs((ynew - y)/ynew);
        if Error < tol
            break
        else y = ynew;
            q = q+1;
        end 
        R_crit(i) = r_crit(abs(ynew));
    end 
end
display(R_crit)
figure;
plot(t,R_crit,'-ro')
legend('R_{crit}')
title('Critical Radii at varying z distances')
ylabel('Radius (cm)')
xlabel('Axial Distance (cm)')
%% 10 Fold Increase

clear all;
clc;

D = 1*10^(-7); %cm^2/s Diffusion Coefficient
V = 0.002; %cm/s Velocity
rc = 0.0005; %cm capillary radius
tm = 5*10^(-5); %cm capillary wall thickness
Ko = 5.75*10^(-5); %cm/s Overall Mass Transfer Rate
Ro1 = 0.1; %umole/cm^3 * s Metabolic Consumption Rate
Co = 5; %umole/cm^3
t = [0.001:0.01:0.101]; %varying z from 0.001 to 0.1 by 0.01

r_crit = @(y) y.*(rc + tm);
Krogh = @(y,z) y.^2.*log(y.^2) - y.^2 + 1 - ((4*D*Co)/(Ro1*(rc+tm)^2)) + 4*(D/(V*rc^2)).*(y.^2-1)*z+(2*D/(rc*Ko)).*(y.^2-1);
dKdy = @(y,z) 2*y+2.*y.*log(y.^2) - 2*y +8*(D/(V*rc^2)).*y.*z + (4*D/(rc*Ko)).*y;

tol = 10^(-7);
maxit=100;
for i=1:11
    z=t(i);
    y=0.001; q=1;
    while q<100
        ynew = y - (Krogh(y,z)/dKdy(y,z));
        Error = abs((ynew - y)/ynew);
        if Error < tol
            break
        else y = ynew;
            q = q+1;
        end 
        R_crit1(i) = r_crit(abs(ynew));
    end 
end
R_crit1
figure;
plot(t,R_crit1,'-bo')
legend('R_{crit}')
title('Critical Radii at varying z distances')
ylabel('Radius (cm)')
xlabel('Axial Distance (cm)')

%% Combined Information

clear all;
clc;

D = 1*10^(-7); %cm^2/s Diffusion Coefficient
V = 0.002; %cm/s Velocity
rc = 0.0005; %cm capillary radius
tm = 5*10^(-5); %cm capillary wall thickness
Ko = 5.75*10^(-5); %cm/s Overall Mass Transfer Rate
Ro = 0.01; %umole/cm^3 * s Metabolic Consumption Rate
Ro1 = 0.1; %umole/cm^3 * s Metabolic Consumption Rate
Co = 5; %umole/cm^3
t = [0.001:0.01:0.101]; %varying z from 0.001 to 0.1 by 0.01

r_crit = @(y) y.*(rc + tm);
Krogh = @(y,z) y.^2.*log(y.^2) - y.^2 + 1 - ((4*D*Co)/(Ro*(rc+tm)^2)) + 4*(D/(V*rc^2)).*(y.^2-1)*z+(2*D/(rc*Ko)).*(y.^2-1);
dKdy = @(y,z) 2.*y+2.*y.*log(y.^2) - 2*y +8*(D/(V*rc^2)).*y.*z + (4*D/(rc*Ko)).*y;

Krogh1 = @(y,z) y.^2.*log(y.^2) - y.^2 + 1 - ((4*D*Co)/(Ro1*(rc+tm)^2)) + 4*(D/(V*rc^2)).*(y.^2-1)*z+(2*D/(rc*Ko)).*(y.^2-1);
dKdy1 = @(y,z) 2.*y+2.*y.*log(y.^2) - 2*y +8*(D/(V*rc^2)).*y.*z + (4*D/(rc*Ko)).*y;

tol = 10^(-7);
maxit=100;
for i=1:11
    z=t(i);
    y=0.001; q=1;
    while q<100
        ynew = y - (Krogh(y,z)/dKdy(y,z));
        Error = abs((ynew - y)/ynew);
        if Error < tol
            break
        else y = ynew;
            q = q+1;
        end 
        R_crit(i) = r_crit(abs(ynew));
    end 
end
display(R_crit)

for i=1:11
    z=t(i);
    y=0.001; q=1;
    while q<100
        ynew = y - (Krogh1(y,z)/dKdy1(y,z));
        Error = abs((ynew - y)/ynew);
        if Error < tol
            break
        else y = ynew;
            q = q+1;
        end 
        R_crit1(i) = r_crit(abs(ynew));
    end 
end
R_crit1
figure;
plot(t,R_crit,'-ro', t, R_crit1,'-bo')
legend('R_{crit} when R = 0.01','R_{crit} when R = 0.1')
title('Critical Radii at varying z distances')
ylabel('Radius (cm)')
xlabel('Axial Distance (cm)')

Z = [0.001:0.01:0.101]';
R_crit_Initial = R_crit';
R_crit_TenFold = R_crit1';

Table = table(Z,R_crit_Initial,R_crit_TenFold)

%% Comparison

%As shown by the graphs and the data table in the sections above, there is
%a clear difference between when Ro = 0.01 and when Ro = 0.1. Ro represents
% the "Metabolic Consumption Rate", meaning how quickly the tissues consume
% the solute brought before them. As a result, we would expect that if the
% tissues consume the solute much more quickly, the solute will not be able
% to travel as far. This hypothesis is confirmed when analyzing the graph
% between the two Ro values. The graph in which Ro = 0.01 has a much higher
% critical radius for each z value than that of the graph in which Ro =
% 0.1. Mathematically, the the difference between the two graphs is about
% 3-fold, wherein Ro=0.01 has a critical radius 3 times as great as the
% critical radius of the graph when Ro=0.1, at laest in the beginning. As z
% increases, the difference between the two graphs lessens, however the
% Ro=0.01 graph always maintains a higher critical radius at every z value
% than the Ro=0.1 graph. Again, this is to be expected because the slower
% tissue consumption rate allows the solute to travel further into the
% tissue, hence a greater critical radius. 
