L = 0.02;
delta_x = 0;
% Sc and St are the coefficients from Q = Sc  +  St * T
N_points_interior = 600;
Thermal_cond = 0.5;
Q_gen = 1000000; % W/m^3
area_of_crossection = 1;
St = -1;
Sc = 1000;
T_init = 100;
Perimeter = 1;
Conv_coefficient = 0;
Surrounding_temp = 200;
Z = Get_slab_temperature(Thermal_cond,Q_gen, T_init, Sc, St, N_points_interior, L, area_of_crossection,Conv_coefficient, Perimeter, Surrounding_temp);
Y = Get_slab_analytical(Q_gen, Thermal_cond,L,N_points_interior,Surrounding_temp, T_init);
plot(Z(:,2),Z(:,1), Color=[0 1 1], LineWidth=5);
xlabel('X');
ylabel('Y');
title('Temperature Variations');
hold("on");
plot(Y(:,2),Y(:,1), Color=[1 0 0], LineWidth=2);
xlabel('X');
ylabel('Y');
title('Temperature Variations Analytical');

function [array_temps_analytical] = Get_slab_analytical(Q_flux, K,L,N_points,T_surr, T_b)
    array_temps_analytical = zeros(N_points+2, 2);
    delta_x = zeros(N_points + 2, 1);

    for i = 2:1:N_points+1
        delta_x(i) = L/N_points;
    end

    delta_x(1) = L/(N_points * 2);
    delta_x(N_points+2) = L/(N_points * 2);
    x(1) = 0;
    for i = 2:1:N_points+2 % Getting all the points from the delta
        x(i) = x(i-1) + delta_x(i-1);
    end

    % array_temps_analytical(1,1) = T_b; 
    for i = 1:1:N_points+2
        array_temps_analytical(i,1) =  T_b + ((T_surr - T_b)/L)*x(i) + (Q_gen*L/(2*K))*x(i) - Q_gen*x(i)*x(i)/(2*K);
        array_temps_analytical(i,2) = x(i);
        % disp(array_temps(i))
    end

    
end

% disp(Z);

function [array_temps] = Get_slab_temperature(Thermal_cond,Q_gen, T_init, Sc, St, N_points, L, area, h, p, T_surr)

% Initializing important parameters for the given problem

x = zeros(N_points+2, 1); % For the unit lengths
a = zeros(N_points+2, 1); % For the right side of the matrix diagonal elements
b = zeros(N_points+2, 1); % For the left side of the matrix diagonal elements
d = zeros(N_points+2, 1); % For the diagonal Elements of the matrix
c = zeros(N_points+2, 1); % For the equation's answers
array_temps = zeros(N_points+2, 2); % For the equation's answers
delta_x = zeros(N_points + 2, 1);

for i = 2:1:N_points+1
    delta_x(i) = L/N_points;
end

delta_x(1) = L/(N_points * 2);
delta_x(N_points+2) = L/(N_points * 2);
% disp(delta_x)

x(1) = 0;
for i = 2:1:N_points+2 % Getting all the points from the delta
    x(i) = x(i-1) + delta_x(i-1);
end


% Initial Boundary Conditions

a(1) = 0;
b(1) = 0;
c(1) = T_init;
d(1) = 1;

m_squared = h*p/(Thermal_cond*area);

% Thermal Conductivity for the remaining points except for the last one
for i = 2:1:N_points+1
    % disp(i);
    d(i) = -(1/(delta_x(i-1)) + 1/delta_x(i+1) + m_squared*(delta_x(i+1) + delta_x(i-1))); % here's where St comes in use
    % disp(d(i));
    a(i) = (1/delta_x(i+1));
    b(i) = (1/delta_x(i-1));
    c(i) = (-m_squared*T_surr- Q_gen)*(delta_x(i+1) + delta_x(i-1)); % here's where Sc comes in use
end

% Final Boundary Conditions

a(N_points + 2) = 0;
b(N_points + 2) = 0;
d(N_points + 2) = 1;
c(N_points + 2) = T_surr;


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
