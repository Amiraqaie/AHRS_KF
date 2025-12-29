function x_dot = state_derivative(x, omega_m)
    phi   = x(1);
    theta = x(2);
    psi   = x(3);
    b     = x(4:6);
    
    % true angular velocity
    omega = omega_m - b;
    p = omega(1);
    q = omega(2);
    r = omega(3);
    
    % Transformation matrix from body rates to Euler rates
    T = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
         0, cos(phi),           -sin(phi);
         0, sin(phi)/cos(theta), cos(phi)/cos(theta)];
    
    % Euler angles derivative
    euler_dot = T * omega;
    
    % Gyro bias derivative (assume random walk = 0)
    b_dot = zeros(3,1);
    
    x_dot = [euler_dot; b_dot];
end
