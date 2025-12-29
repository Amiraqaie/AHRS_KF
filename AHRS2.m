clc;
clear;
close all;
jacobian;  % Make sure your functions Fk, Lk, state_derivative exist

%% Load dataset
Dataset = readtable('Dataset2.csv');
Ts = 0.0099;

%% Kalman Filter parameters
Q = 5e-1 * eye(3);      
R = diag([2.5e-3, 2.5e-3, 1e-2]);                   
P = eye(6);            
xplus = zeros(6,1);

C = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0];

N = height(Dataset);

% Preallocate arrays
states = zeros(N,6);
P_store = zeros(6,6,N);
meas = zeros(N,3);

%% Kalman Filter loop
for k = 1:N
    % --- Read measurements ---
    ax = Dataset.Acceleration_X(k);  % already in m/s²
    ay = Dataset.Acceleration_Y(k);
    az = Dataset.Acceleration_Z(k);

    omega = deg2rad([Dataset.Gyro_X_dps_(k); Dataset.Gyro_Y_dps_(k); Dataset.Gyro_Z_dps_(k)]); % deg/s → rad/s

    mx = Dataset.Magnetometer_X(k); % already in µT
    my = Dataset.Magnetometer_Y(k);
    mz = Dataset.Magnetometer_Z(k);

    % --- Euler angles from sensors ---
    phi_ = atan2(ay, az);
    theta_ = atan2(-ax, sqrt(ay^2 + az^2));
    psi_ = atan2(mz*sin(phi_) - my*cos(phi_), ...
                 mx*cos(theta_) + my*sin(theta_)*sin(phi_) + mz*sin(theta_)*cos(phi_));

    Zt = [phi_; theta_; psi_];
    meas(k,:) = Zt;

    if k == 1
        xplus = [Zt; 0; 0; 0];
        states(k,:) = xplus;
        P_store(:,:,k) = P;
        continue
    end

    % --- Predict ---
    F_k = Fk(xplus, omega, Ts);
    L_k = Lk(xplus, omega, Ts);
    x_ = xplus + state_derivative(xplus, omega) * Ts;
    P = F_k*P*F_k' + L_k*Q*L_k';

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