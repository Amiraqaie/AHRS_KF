clc; clear; close all;

Dataset = readtable('Dataset1.xlsx');
Ts = 0.0099;

beta = 0.05;
N = height(Dataset);

q = [1 0 0 0]';

euler_est = zeros(N,3);

for k = 1:N

    ax_g = Dataset.Acceleration_X(k) / 1000;
    ay_g = Dataset.Acceleration_Y(k) / 1000;
    az_g = Dataset.Acceleration_Z(k) / 1000;

    ax = ax_g * 9.81;
    ay = ay_g * 9.81;
    az = az_g * 9.81;

    acc = [
        ax
        ay
        az
    ];

    gyro = [Dataset.Gyro_X_mrps_(k); Dataset.Gyro_Y_mrps_(k); Dataset.Gyro_Z_mrps_(k)] / 1000;


    mx = Dataset.Magnetometer_X(k) * 0.1;
    my = Dataset.Magnetometer_Y(k) * 0.1;
    mz = Dataset.Magnetometer_Z(k) * 0.1;

    mag = [
        mx
        my
        mz
    ];

    q = madgwick_update(q, gyro, acc, mag, beta, Ts);

    euler_est(k,:) = quat2euler(q)';
end

figure;
plot(euler_est);
legend('\phi','\theta','\psi');
title('Madgwick AHRS');
