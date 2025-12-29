clc;
clear;
close all;

%% Load dataset
Dataset = readtable('Dataset1.xlsx');
Ts = 0.0203;  % Sampling time

%% Kalman Filter Parameters
Q = 5e-1 * eye(3);     
R = diag([2.5e-3, 2.5e-3, 1e-2]);                  
P = 1e3*eye(6);            
xplus = zeros(6,1);           

C = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0];            % Measurement matrix

N = height(Dataset);

% Preallocate arrays
states = zeros(N, 6);
P_store = zeros(6,6,N);
meas = zeros(N,3);

%% Kalman Filter Loop
for k = 1:N
    % --- Read sensor measurements ---
    ax_g = Dataset.Acceleration_X(k) / 1000;
    ay_g = Dataset.Acceleration_Y(k) / 1000;
    az_g = Dataset.Acceleration_Z(k) / 1000;

    ax = ax_g * 9.81;
    ay = ay_g * 9.81;
    az = az_g * 9.81;

    omega = [Dataset.Gyro_X_mrps_(k); Dataset.Gyro_Y_mrps_(k); Dataset.Gyro_Z_mrps_(k)] / 1000;

    mx = Dataset.Magnetometer_X(k) * 0.1;
    my = Dataset.Magnetometer_Y(k) * 0.1;
    mz = Dataset.Magnetometer_Z(k) * 0.1;

    % --- Calculate Euler angles from sensors ---
    phi_ = atan2(ay, az);
    theta_ = atan2(-ax, sqrt(ay^2 + az^2));
    psi_ = atan2(mz*sin(phi_) - my*cos(phi_), ...
                 mx*cos(theta_) + my*sin(theta_)*sin(phi_) + mz*sin(theta_)*cos(phi_));

    Zt = [phi_; theta_; psi_];
    meas(k,:) = Zt;

    if k == 1
        xplus = [Zt; 0; 0; 0];  % Initialize state with angles + zeros for angular velocity
        states(k,:) = xplus;
        P_store(:,:,k) = P;
        continue
    end

    % --- Compute Jacobians ---
    F_k = Fk(xplus, omega, Ts);  % User-defined function
    L_k = Lk(xplus, omega, Ts);  % User-defined function

    % --- Predict ---
    P = F_k*P*F_k' + L_k*Q*L_k';
    x_ = xplus + state_derivative(xplus, omega) * Ts;  % User-defined function

    % --- Update ---
    S = C*P*C' + R;
    K = (P*C') / S;

    innov = Zt - C*x_;
    for i = 1:3
        while innov(i) > pi
            innov(i) = innov(i) - 2*pi;
        end
        while innov(i) < -pi
            innov(i) = innov(i) + 2*pi;
        end
    end
    xplus = x_ + K*innov;

    P = (eye(6) - K*C) * P;

    % --- Store results ---
    states(k,:) = xplus;
    P_store(:,:,k) = P;
end
result;
