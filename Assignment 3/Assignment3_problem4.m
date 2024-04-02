
% Lengths of the two one-dimensional sections
La = 0.05;

% Left face of the wall temperature in *C
T_init = 30;
T_right = 600;

% For simplicity let's keep the grid sizes to be constant
delta_x = 0;


N_points_interior = 100;

Sc = 0 ;
Sp = 0 ;

Z = Get_slab_temperature(T_init,T_right,La, N_points_interior, Sc, Sp);

plot(Z(:,2),Z(:,1), Color=[0 0 1], LineWidth=5);
xlabel('X');
ylabel('Y');
title('Temperature Variations Numerical for Conjunction of both rods');

function [array_temps] = Get_slab_temperature(T_init,T_right,La, N_points,Sc,Sp)

% Initializing important parameters for the given problem

x = zeros(N_points+2, 1); % For the unit lengths
a = zeros(N_points+2, 1); % For the right side of the matrix diagonal elements
b = zeros(N_points+2, 1); % For the left side of the matrix diagonal elements
d = zeros(N_points+2, 1); % For the diagonal Elements of the matrix
c = zeros(N_points+2, 1); % For the equation's answers
array_temps = zeros(N_points+2, 2); % For the equation's answers
array_temps_guesses = zeros(N_points+2, 2);
delta_x = zeros(N_points + 1, 1);
% delta_x_b = zeros(N_points + 1, 1);

for i = 2:1:N_points+1
    delta_x(i) = (La)/N_points;
    % delta_x_b(i) = (Lb)/N_points;
end

delta_x(1) = (La)/(N_points * 2);
delta_x(N_points + 2) = La/(N_points*2);

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

a(1) = 0;
b(1) = 0;
d(1) = 1;
c(1) = T_init;
m_squared = 0;
    
% Final Boundary Conditions

a(N_points + 2) = 0;
b(N_points + 2) = 0;
d(N_points + 2) = 1;
c(N_points + 2) = T_right;
    

error_acceptable = 0.001;
error = T_right;
    while error > error_acceptable
    % Thermal Conductivity for the remaining points except for the last one
        for i = 2:1:N_points+1
            % disp(i);
        
            a(i) = -((1 + 0.008*array_temps_guesses(i-1,1))/delta_x(i+1));
            b(i) = -((1 + 0.008*array_temps_guesses(i+1,1))/delta_x(i-1));    
            d(i) = -( a(i) + b(i));
            % disp(d(i));
            c(i) =  0; % here's where Sc comes in use
        end
        
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
        end
        error = sqrt(mean(array_temps(:,1) - array_temps_guesses(:,1)).^2);
        disp(error);
        array_temps_guesses = array_temps;
    end

    disp("Solution Found !!!\n");
    disp(array_temps);

end