% Create filled contour plot from csv-file with the solution

rho = readmatrix("output/rho_0.csv")';
u   = readmatrix("output/u_0.csv")';
v   = readmatrix("output/v_0.csv")';
w   = readmatrix("output/w_0.csv")';
p   = readmatrix("output/p_0.csv")';
T   = readmatrix("output/T_0.csv")';

x = linspace(0,5,41);
z = linspace(0,1,41);

norm_rho   = readmatrix("output/norm_rho.dat");
norm_rho_u = readmatrix("output/norm_rho_u.dat");
norm_rho_v = readmatrix("output/norm_rho_v.dat");
norm_rho_w = readmatrix("output/norm_rho_w.dat");
norm_E     = readmatrix("output/norm_E.dat");

figure(1);
[M, myPlot] = contourf(x, z, u, 30);
%myPlot.LineColor = 'None';
colorbar;
title('Velocity component u');
xlabel('x');
ylabel('z');

figure(2);
hold on
plot(u(:,1),   z, 'DisplayName', 'Inlet (analytic)');
plot(u(:,end), z, 'DisplayName', 'Outlet');
legend;
title('Velocity component u');
xlabel('u');
ylabel('z');

figure(3);
[M, myPlot] = contourf(x, z, p, 30);
%myPlot.LineColor = 'None';
colorbar;
title('Pressure, p');
xlabel('x');
ylabel('z');

figure(4);
plot(x, p(21,:));
title('Pressure, p');
xlabel('x');
ylabel('p');

figure(5);
plot(x, rho(21,:));
title('Density, \rho');
xlabel('x');
ylabel('\rho');

figure(6);
[M, myPlot] = contourf(x, z, w, 30);
%myPlot.LineColor = 'None';
colorbar;
title('Velocity component w');
xlabel('x');
ylabel('z');

figure(7);
[M, myPlot] = contourf(x, z, T, 30);
%myPlot.LineColor = 'None';
colorbar;
title('Temperature, T');
xlabel('x');
ylabel('z');

set(0,'DefaultLineLineWidth',2)

figure(8);
plot(norm_rho(2:end) / norm_rho(2));
title('Norm history, \rho');
xlabel('time level, n');
ylabel('||\rho||');
set(gca, 'YScale', 'log')

figure(9);
plot(norm_rho_u(2:end) / norm_rho_u(2));
title('Norm history, \rhou');
xlabel('time level, n');
ylabel('||\rhou||');
set(gca, 'YScale', 'log')

figure(10);
plot(norm_rho_v(2:end) / norm_rho_v(2));
title('Norm history, \rhov');
xlabel('time level, n');
ylabel('||\rhov||');
set(gca, 'YScale', 'log')

figure(11);
plot(norm_rho_w(2:end) / norm_rho_w(2));
title('Norm history, \rhow');
xlabel('time level, n');
ylabel('||\rhow||');
set(gca, 'YScale', 'log')

figure(12);
plot(norm_E(2:end) / norm_E(2));
title('Norm history, \rhoE');
xlabel('time level, n');
ylabel('||\rhoE||');
set(gca, 'YScale', 'log')

fprintf("max(v) = %1.2f \n", max(max(v)));
fprintf("min(v) = %1.2f \n", min(min(v)));
fprintf( "Density ratio: %1.4f \n", (rho(21,1)+1)/(rho(21,end)+1));
fprintf("Velocity ratio: %1.4f \n", (u(21,end))/(u(21,1)));




