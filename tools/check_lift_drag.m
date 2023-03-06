clear;
close all;

NI=1600;
NJ=802;
NK=3;
nNodes = NI*NJ*NK;

L_x = 40;
L_y = 20;
x_c = L_x/4;
y_c = L_y/2;
D = 1;
r_c = D/2;
Ma = 0.25;
u_ref = Ma;
rho_0 = 1;
gamma = 1.4;
Re = 40/Ma;
ScPlusOne = 1 + 110.4/300;

formatString = "%s";
formatFloat = '%f';
formatInt = '%i';

filePath = "/media/frederk/Dump drive/Simulations/Cylinder extrap Re40 Ma0.25 1600x802x3/output/out.vtk.9";

fileID = fopen(filePath, 'r');
fscanf(fileID, formatString, 34);
flag_1D=fscanf(fileID, formatInt, nNodes); % NodeFlag
fscanf(fileID, formatString, 6);
fscanf(fileID, formatFloat, nNodes); % density
fscanf(fileID, formatString, 6);
p_1D=fscanf(fileID, formatFloat, nNodes); % pressure
fscanf(fileID, formatString, 6);
T_1D=fscanf(fileID, formatFloat, nNodes); % Temperature
fscanf(fileID, formatString, 3);
vel_1D=fscanf(fileID, formatFloat, nNodes*3); % Velocity
fclose(fileID);

flag = reshape(flag_1D, [NI,NJ,NK]);
flag = reshape(flag(:,2:end-1,2), [NI,NJ-2]);
p = reshape(p_1D, [NI,NJ,NK]);
p = reshape(p(:,2:end-1,2), [NI,NJ-2]);
T = reshape(T_1D, [NI,NJ,NK]);
T = reshape(T(:,2:end-1,2), [NI,NJ-2]);
velocity = reshape(vel_1D, [3,NI,NJ,NK]);
u = reshape(velocity(1,:,2:end-1,2), [NI, NJ-2]);
v = reshape(velocity(2,:,2:end-1,2), [NI, NJ-2]);
w = reshape(velocity(3,:,2:end-1,2), [NI, NJ-2]);
x = linspace(0,L_x,NI);
y = linspace(0,L_y,NJ-2);
dx = L_x/(NI-1);
dy = L_y/(NJ-3);

BI_points = [];
angles = [];
normal_vecs = [];
traction_vecs = [];
ghostNodes = find(flag==3)';
for index1D = ghostNodes
    j = floor( (index1D-1) / NI) + 1;
    i = rem(index1D-1, NI) + 1;
    if flag(i+1,j) == 0
        y_BI = y(j);
        x_BI = x_c + sqrt(r_c^2 - (y_BI-y_c)^2);
        BI = [x_BI; y_BI];
        BI_points(1:2,end+1) = BI;
        angle_BI = atan2(y_BI-y_c, x_BI-x_c);
        if angle_BI < 0
            angle_BI = angle_BI + 2*pi;
        end
        angles(end+1) = angle_BI;
        p_BI = p(i,j)+(p(i+1,j)-p(i,j))/dx*(x_BI-x(i));
        T_BI = T(i,j)+(T(i+1,j)-T(i,j))/dx*(x_BI-x(i));
        mu_BI = (1+T_BI)^1.5 * ScPlusOne / ( Re*( T_BI + ScPlusOne ) );
        dudx_1 = (u(i+2,j)-u(i,j))/(2*dx);
        dudx_2 = (u(i+3,j)-u(i+1,j))/(2*dx);
        dudx_BI = dudx_1+(dudx_1-dudx_2)/-dx*(x_BI-x(i+1));
        dvdx_1 = (v(i+2,j)-v(i,j))/(2*dx);
        dvdx_2 = (v(i+3,j)-v(i+1,j))/(2*dx);
        dvdx_BI = dvdx_1+(dvdx_1-dvdx_2)/-dx*(x_BI-x(i+1));
        dudy_1 = (u(i+1,j+1)-u(i+1,j-1))/(2*dy);
        dudy_2 = (u(i+2,j+1)-u(i+2,j-1))/(2*dy);
        dudy_BI = dudy_1+(dudy_1-dudy_2)/-dx*(x_BI-x(i+1));
        dvdy_1 = (v(i+1,j+1)-v(i+1,j-1))/(2*dy);
        dvdy_2 = (v(i+2,j+1)-v(i+2,j-1))/(2*dy);
        dvdy_BI = dvdy_1+(dvdy_1-dvdy_2)/-dx*(x_BI-x(i+1));
		deformation = [dudx_BI, (dudy_BI+dvdx_BI)/2; (dudy_BI+dvdx_BI)/2, dvdy_BI ];
		viscousStress = 2*mu_BI*deformation - 2./3.*mu_BI*(dudx_BI+dvdy_BI)*eye(2);
		totalStress = viscousStress - p_BI*eye(2);
		normal = [x_BI-x_c; y_BI-y_c];
		normal = normal ./ norm(normal);
        normal_vecs(1:2, end+1) = normal;
        traction_vecs(1:2, end+1) = totalStress * normal;
    end
    if flag(i-1,j) == 0
        y_BI = y(j);
        x_BI = x_c - sqrt(r_c^2 - (y_BI-y_c)^2);
        BI = [x_BI; y_BI];
        BI_points(1:2,end+1) = BI;
        angle_BI = atan2(y_BI-y_c, x_BI-x_c);
        if angle_BI < 0
            angle_BI = angle_BI + 2*pi;
        end
        angles(end+1) = angle_BI;
        p_BI = p(i,j)+(p(i-1,j)-p(i,j))/-dx*(x_BI-x(i));
        T_BI = T(i,j)+(T(i-1,j)-T(i,j))/-dx*(x_BI-x(i));
        mu_BI = (1+T_BI)^1.5 * ScPlusOne / ( Re*( T_BI + ScPlusOne ) );
        dudx_1 = (u(i,j)-u(i-2,j))/(2*dx);
        dudx_2 = (u(i-1,j)-u(i-3,j))/(2*dx);
        dudx_BI = dudx_1+(dudx_1-dudx_2)/dx*(x_BI-x(i-1));
        dvdx_1 = (v(i,j)-v(i-2,j))/(2*dx);
        dvdx_2 = (v(i-1,j)-v(i-3,j))/(2*dx);
        dvdx_BI = dvdx_1+(dvdx_1-dvdx_2)/dx*(x_BI-x(i-1));
        dudy_1 = (u(i-1,j+1)-u(i-1,j-1))/(2*dy);
        dudy_2 = (u(i-2,j+1)-u(i-2,j-1))/(2*dy);
        dudy_BI = dudy_1+(dudy_1-dudy_2)/dx*(x_BI-x(i-1));
        dvdy_1 = (v(i-1,j+1)-v(i-1,j-1))/(2*dy);
        dvdy_2 = (v(i-2,j+1)-v(i-2,j-1))/(2*dy);
        dvdy_BI = dvdy_1+(dvdy_1-dvdy_2)/dx*(x_BI-x(i-1));
		deformation = [dudx_BI, (dudy_BI+dvdx_BI)/2; (dudy_BI+dvdx_BI)/2, dvdy_BI ];
		viscousStress = 2*mu_BI*deformation - 2./3.*mu_BI*(dudx_BI+dvdy_BI)*eye(2);
		totalStress = viscousStress - p_BI*eye(2);
		normal = [x_BI-x_c; y_BI-y_c];
		normal = normal ./ norm(normal);
        normal_vecs(1:2, end+1) = normal;
        traction_vecs(1:2, end+1) = totalStress * normal;
    end
    if flag(i,j+1) == 0
        x_BI = x(i);
        y_BI = y_c + sqrt(r_c^2 - (x_BI-x_c)^2);
        BI = [x_BI; y_BI];
        BI_points(1:2,end+1) = BI;
        angle_BI = atan2(y_BI-y_c, x_BI-x_c);
        if angle_BI < 0
            angle_BI = angle_BI + 2*pi;
        end
        angles(end+1) = angle_BI;
        p_BI = p(i,j)+(p(i,j+1)-p(i,j))/dy*(y_BI-y(j));
        T_BI = T(i,j)+(T(i,j+1)-T(i,j))/dy*(y_BI-y(j));
        mu_BI = (1+T_BI)^1.5 * ScPlusOne / ( Re*( T_BI + ScPlusOne ) );
        dudx_1 = (u(i+1,j+1)-u(i-1,j+1))/(2*dx);
        dudx_2 = (u(i+1,j+2)-u(i-1,j+2))/(2*dx);
        dudx_BI = dudx_1+(dudx_1-dudx_2)/-dy*(y_BI-y(j+1));
        dvdx_1 = (v(i+1,j+1)-v(i-1,j+1))/(2*dx);
        dvdx_2 = (v(i+1,j+2)-v(i-1,j+2))/(2*dx);
        dvdx_BI = dvdx_1+(dvdx_1-dvdx_2)/-dy*(y_BI-y(j+1));
        dudy_1 = (u(i,j+2)-u(i,j))/(2*dy);
        dudy_2 = (u(i,j+3)-u(i,j+1))/(2*dy);
        dudy_BI = dudy_1+(dudy_1-dudy_2)/-dy*(y_BI-y(j+1));
        dvdy_1 = (v(i,j+2)-v(i,j))/(2*dy);
        dvdy_2 = (v(i,j+3)-v(i,j+1))/(2*dy);
        dvdy_BI = dvdy_1+(dvdy_1-dvdy_2)/-dy*(y_BI-y(j+1));
		deformation = [dudx_BI, (dudy_BI+dvdx_BI)/2; (dudy_BI+dvdx_BI)/2, dvdy_BI ];
		viscousStress = 2*mu_BI*deformation - 2./3.*mu_BI*(dudx_BI+dvdy_BI)*eye(2);
		totalStress = viscousStress - p_BI*eye(2);
		normal = [x_BI-x_c; y_BI-y_c];
		normal = normal ./ norm(normal);
        normal_vecs(1:2, end+1) = normal;
        traction_vecs(1:2, end+1) = totalStress * normal;
    end
    if flag(i,j-1) == 0
        x_BI = x(i);
        y_BI = y_c - sqrt(r_c^2 - (x_BI-x_c)^2);
        BI = [x_BI; y_BI];
        BI_points(1:2,end+1) = BI;
        angle_BI = atan2(y_BI-y_c, x_BI-x_c);
        if angle_BI < 0
            angle_BI = angle_BI + 2*pi;
        end
        angles(end+1) = angle_BI;
        p_BI = p(i,j)+(p(i,j-1)-p(i,j))/-dy*(y_BI-y(j));
        T_BI = T(i,j)+(T(i,j-1)-T(i,j))/-dy*(y_BI-y(j));
        mu_BI = (1+T_BI)^1.5 * ScPlusOne / ( Re*( T_BI + ScPlusOne ) );
        dudx_1 = (u(i+1,j-1)-u(i-1,j-1))/(2*dx);
        dudx_2 = (u(i+1,j-2)-u(i-1,j-2))/(2*dx);
        dudx_BI = dudx_1+(dudx_1-dudx_2)/dy*(y_BI-y(j-1));
        dvdx_1 = (v(i+1,j-1)-v(i-1,j-1))/(2*dx);
        dvdx_2 = (v(i+1,j-2)-v(i-1,j-2))/(2*dx);
        dvdx_BI = dvdx_1+(dvdx_1-dvdx_2)/dy*(y_BI-y(j-1));
        dudy_1 = (u(i,j)-u(i,j-2))/(2*dy);
        dudy_2 = (u(i,j-1)-u(i,j-3))/(2*dy);
        dudy_BI = dudy_1+(dudy_1-dudy_2)/dy*(y_BI-y(j-1));
        dvdy_1 = (v(i,j)-v(i,j-2))/(2*dy);
        dvdy_2 = (v(i,j-1)-v(i,j-3))/(2*dy);
        dvdy_BI = dvdy_1+(dvdy_1-dvdy_2)/dy*(y_BI-y(j-1));
		deformation = [dudx_BI, (dudy_BI+dvdx_BI)/2; (dudy_BI+dvdx_BI)/2, dvdy_BI ];
		viscousStress = 2*mu_BI*deformation - 2./3.*mu_BI*(dudx_BI+dvdy_BI)*eye(2);
		totalStress = viscousStress - p_BI*eye(2);
		normal = [x_BI-x_c; y_BI-y_c];
		normal = normal ./ norm(normal);
        normal_vecs(1:2, end+1) = normal;
        traction_vecs(1:2, end+1) = totalStress * normal;
    end
end
[angles, sortIndex] = sort(angles);
BI_points = BI_points(:, sortIndex);
normal_vecs = normal_vecs(:, sortIndex);
traction_vecs = traction_vecs(:, sortIndex);

prev_BI = BI_points(:, end);
prev_traction = traction_vecs(:, end);
force = [0; 0];

for k=1:length(angles)
    current_BI = BI_points(:, k);
    current_surface = norm( current_BI - prev_BI );
    current_traction = traction_vecs(:, k);
    force_contrib = (current_traction + prev_traction)/2 * current_surface;
    force = force + force_contrib;
    prev_BI = current_BI;
    prev_traction = current_traction;
end

C_D = force(1) * 2 / ( rho_0 * u_ref^2 * D );
C_L = force(2) * 2 / ( rho_0 * u_ref^2 * D );

figure(1);
plot(BI_points(1,:), BI_points(2,:));
hold on;
quiver(BI_points(1,:), BI_points(2,:), traction_vecs(1,:), traction_vecs(2,:));
axis equal;


