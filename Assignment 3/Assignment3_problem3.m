
% Lengths of the two one-dimensional sections
La = 0.4;
Lb = 0.2;

% Thermal conductivities of the two one-dimentional sections
Ka = 10 ; % W/m-K
Kb = 0.1; % W/m-K

% Left face of the wall temperature in *C
T_init = 800;
T_right = 20;

% Heat transfer coefficient over the left surface
h_left = 25; % W/m^2-K
h_right = 10; % W/m^2-K

% For simplicity let's keep the grid sizes to be constant
delta_x = 0;

N_points_interior = 100;

Sc = 0 ;
Sp = 0 ;

Z = Get_slab_temperature(Ka,Kb,T_init,T_right,La,Lb,h_left,h_right, N_points_interior, Sc, Sp);
Y = Get_slab_analytical(Ka,Kb,T_init,T_right,La,Lb,h_left,h_right, N_points_interior);


plot(Z(:,2),Z(:,1), Color=[0 0 1], LineWidth=5);
xlabel('X');
ylabel('Y');
title('Temperature Variations Numerical for Conjunction of both rods');


hold("on");
plot(Y(:,2),Y(:,1), Color=[0 1 0], LineWidth=2);
xlabel('X');
ylabel('Y');
title('Temperature Variations Analytical');

function [array_temps_analytical] = Get_slab_analytical(Ka,Kb,T_init,T_right,La,Lb,h_left,h_right, N_points)
    array_temps_analytical = zeros(N_points+2, 2);
    delta_x = zeros(N_points + 2, 1);

    R_equivalent = ((La/Ka) + (Lb/Kb) + (1/h_left) + 1/(h_right));
    disp(R_equivalent);

    L = La + Lb;
    for i = 2:1:N_points+1
        delta_x(i) = L/N_points;
    end

    delta_x(1) = L/(N_points * 2);
    delta_x(N_points+2) = L/(N_points * 2);
    x(1) = 0;
    for i = 2:1:N_points+2 % Getting all the points from the delta
        x(i) = x(i-1) + delta_x(i-1);
    end
    disp(x);

    array_temps_analytical(1,1) = ((T_right - T_init)/(h_left*R_equivalent)) + T_init;
    array_temps_analytical(N_points+2,1) = -((T_right - T_init)/(h_right*R_equivalent)) + T_right;
    array_temps_analytical(1,2) = x(1);
    array_temps_analytical(N_points+2,2) = x(N_points+2);

    for i = 2:1:N_points+1
        if x(i) < La
            array_temps_analytical(i,1) = ((T_right - T_init)*(delta_x(i))/(Ka*R_equivalent)) + array_temps_analytical(i-1,1);
        else
            array_temps_analytical(i,1) = ((T_right - T_init)*(delta_x(i))/(Kb*R_equivalent)) + array_temps_analytical(i-1,1);
        end
        array_temps_analytical(i,2) = x(i);
    end

        % disp(array_temps_analytical);
    
end
% disp(Z);

function [array_temps] = Get_slab_temperature(Ka,Kb,T_init,T_right,La,Lb,h_left,h_right,N_points,Sc,Sp)

% Initializing important parameters for the given problem

x = zeros(N_points+2, 1); % For the unit lengths
a = zeros(N_points+2, 1); % For the right side of the matrix diagonal elements
b = zeros(N_points+2, 1); % For the left side of the matrix diagonal elements
d = zeros(N_points+2, 1); % For the diagonal Elements of the matrix
c = zeros(N_points+2, 1); % For the equation's answers
array_temps = zeros(N_points+2, 2); % For the equation's answers
delta_x = zeros(N_points + 1, 1);
% delta_x_b = zeros(N_points + 1, 1);

for i = 2:1:N_points+1
    delta_x(i) = (La+Lb)/N_points;
    % delta_x_b(i) = (Lb)/N_points;
end

delta_x(1) = (La+Lb)/(N_points * 2);
delta_x(N_points + 2) = La/(N_points*2) + Lb/(N_points*2);

% delta_x_b(1) = Lb/(N_points);
% delta_x_b(N_points + 1) = Lb/(N_points*2);


% Making the whole length for the object
x(1) = 0;
for i = 2:1:N_points+2 % Getting all the points from the delta
    x(i) = x(i-1) + delta_x(i-1);
end

% for i = 2:1:N_points+2
%     x(N_points+i) = x(N_points + i - 1) + delta_x_b(i-1);
% end


% Initial Boundary Conditions

a(1) = -Ka/delta_x(1);
b(1) = 0;
d(1) = -a(1) + h_left;
c(1) = h_left*T_init;
m_squared = 0;

% Thermal Conductivity for the remaining points except for the last one
for i = 2:1:N_points+1
    % disp(i);
    
    if x(i+1) < La - 0.02
        d(i) = (Ka/(delta_x(i-1)) + Ka/delta_x(i+1) );
        a(i) = -(Ka/delta_x(i+1));
        b(i) = -(Ka/delta_x(i-1));
    elseif (La - 0.002 <= x(i+1)) && (x(i+1)<= La + 0.002)
        d(i) = (Ka/(delta_x(i-1)) + Kb/delta_x(i+1) );
        a(i) = -(Kb/delta_x(i+1));
        b(i) = -(Ka/delta_x(i-1));
    else
        d(i) = (Kb/(delta_x(i-1)) + Kb/delta_x(i+1) );
        a(i) = -(Kb/delta_x(i+1));
        b(i) = -(Kb/delta_x(i-1));
    end
    
    % disp(d(i));
    c(i) =  0; % here's where Sc comes in use
end

% Final Boundary Conditions

a(N_points + 2) = 0;
b(N_points + 2) = -Kb/delta_x(N_points+2);
d(N_points + 2) = -b(N_points + 2) + h_right;
c(N_points + 2) = h_right*T_right;


for i = 2:1:N_points+2
    d(i) = d(i) - b(i)*a(i-1)/d(i-1);
    c(i) = c(i) - b(i)*c(i-1)/d(i-1);
    
end

% for i = 1:1:N_points+2
% 
%     disp(d(i));
% 
% end

% Final equations for back propogation
array_temps(N_points+2,1) = c(N_points+2)/d(N_points+2);
array_temps(N_points+2,2) = x(N_points+2);
% disp(array_temps(N_points+2))

for i = N_points+1:-1:1
    array_temps(i,1) = (c(i) - a(i)*array_temps(i+1,1))/d(i);
    array_temps(i,2) = x(i);
    disp(array_temps(i));
end


end