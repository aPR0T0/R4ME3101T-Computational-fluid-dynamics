L = 0.5;
K = 1;
delta_x = 0;
% Sc and St are the coefficients from Q = Sc  +  St * T
N_points_interior = 600;
Thermal_cond = 1;
Q_gen = 1000; % W/m^2
area_of_crossection = 1;
St = -1;
Sc = 1000;
T_init = 100;
Z = Get_slab_temperature(Thermal_cond,Q_gen, T_init, Sc, St, N_points_interior, L, area_of_crossection);

plot(Z(:,2),Z(:,1), 'ro-');
xlabel('X');
ylabel('Y');
title('Temperature Variations');


% disp(Z);

function [array_temps] = Get_slab_temperature(Thermal_cond,Q_gen, T_init, Sc, St, N_points, L, area)

% Initializing important parameters for the given problem

x = zeros(N_points+2, 1); % For the unit lengths
a = zeros(N_points+2, 1); % For the right side of the matrix diagonal elements
b = zeros(N_points+2, 1); % For the left side of the matrix diagonal elements
d = zeros(N_points+2, 1); % For the diagonal Elements of the matrix
c = zeros(N_points+2, 1); % For the equation's answers
array_temps = zeros(N_points+2, 2); % For the equation's answers

delta_x = L/(N_points+1);
% disp(delta_x);
x(1) = 0;
for i = 2:1:N_points+2 % Getting all the points from the delta
    x(i) = x(i-1) + delta_x;
end

% Initial Boundary Conditions

a(1) = 0;
b(1) = 0;
c(1) = T_init;
d(1) = 1;

% Thermal Conductivity for the remaining points except for the last one
for i = 2:1:N_points+1
    % disp(i);
    d(i) = -(2-St*((delta_x*delta_x)/Thermal_cond)); % here's where St comes in use
    % disp(d(i));
    a(i) = 1;
    b(i) = 1;
    c(i) = -Sc*(delta_x*delta_x)/Thermal_cond; % here's where Sc comes in use
end

% Final Boundary Conditions

a(N_points + 2) = 0;
b(N_points + 2) = -Thermal_cond*area/delta_x;
d(N_points + 2) = -b(N_points+2);
c(N_points + 2) = Q_gen;


for i = 2:1:N_points+2
    d(i) = d(i) - b(i)*a(i-1)/d(i-1);
    c(i) = c(i) - b(i)*c(i-1)/d(i-1);
    
end

for i = 1:1:N_points+2

    disp(d(i));
    
end

% Final equations for back propogation
array_temps(N_points+2,1) = c(N_points+2)/d(N_points+2);
array_temps(N_points+2,2) = x(N_points+2);
% disp(array_temps(N_points+2))

for i = N_points+1:-1:1
    array_temps(i,1) = (c(i) - a(i)*array_temps(i+1,1))/d(i);
    array_temps(i,2) = x(i);
    % disp(array_temps(i))
end


end

