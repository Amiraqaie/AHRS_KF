clc; clear; close all;

Dataset = readtable('Dataset2.csv');
Ts = 0.0099;

beta = 0.05;     % Typical value
N = height(Dataset);

q = [1 0 0 0]';  % Initial quaternion

euler_est = zeros(N,3);

for k = 1:N
    gyro = deg2rad([
        Dataset.Gyro_X_dps_(k)
        Dataset.Gyro_Y_dps_(k)
        Dataset.Gyro_Z_dps_(k)
    ]);

    acc = [
        Dataset.Acceleration_X(k)
        Dataset.Acceleration_Y(k)
        Dataset.Acceleration_Z(k)
    ];

    mag = [
        Dataset.Magnetometer_X(k)
        Dataset.Magnetometer_Y(k)
        Dataset.Magnetometer_Z(k)
    ];

    q = madgwick_update(q, gyro, acc, mag, beta, Ts);

    euler_est(k,:) = quat2euler(q)';
end

figure;
plot(euler_est);
legend('\phi','\theta','\psi');
title('Madgwick AHRS');
