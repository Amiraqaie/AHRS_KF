%% Extract states for plotting
phi   = states(:,1);
theta = states(:,2);
psi   = states(:,3);
wx    = states(:,4);
wy    = states(:,5);
wz    = states(:,6);

time = Dataset.Time_s_;

%% Plot estimated states
figure;
subplot(2,1,1)
plot(time, phi, 'r', 'LineWidth', 1.2); hold on;
plot(time, theta, 'g', 'LineWidth', 1.2);
plot(time, psi, 'b', 'LineWidth', 1.2);
xlabel('Time [s]')
ylabel('Angles [rad]')
title('Estimated Euler Angles')
legend('\phi','\theta','\psi')
grid on

subplot(2,1,2)
plot(time, wx, 'r', 'LineWidth', 1.2); hold on;
plot(time, wy, 'g', 'LineWidth', 1.2);
plot(time, wz, 'b', 'LineWidth', 1.2);
xlabel('Time [s]')
ylabel('Angular velocities [rad/s]')
title('Estimated Angular Velocities')
legend('\omega_x','\omega_y','\omega_z')
grid on

%% Plot state variances
var_phi   = squeeze(P_store(1,1,:));
var_theta = squeeze(P_store(2,2,:));
var_psi   = squeeze(P_store(3,3,:));
var_wx    = squeeze(P_store(4,4,:));
var_wy    = squeeze(P_store(5,5,:));
var_wz    = squeeze(P_store(6,6,:));

figure;
subplot(2,1,1)
plot(time, var_phi, 'r', 'LineWidth', 1.2); hold on;
plot(time, var_theta, 'g', 'LineWidth', 1.2);
plot(time, var_psi, 'b', 'LineWidth', 1.2);
xlabel('Time [s]')
ylabel('Variance')
title('Variance of Euler Angles')
legend('var(\phi)','var(\theta)','var(\psi)')
grid on
% Set y-axis limits
ylim([0, 5*mean([mean(var_phi), mean(var_theta), mean(var_psi)])]);

subplot(2,1,2)
plot(time, var_wx, 'r', 'LineWidth', 1.2); hold on;
plot(time, var_wy, 'g', 'LineWidth', 1.2);
plot(time, var_wz, 'b', 'LineWidth', 1.2);
xlabel('Time [s]')
ylabel('Variance')
title('Variance of Angular Velocities')
legend('var(\omega_x)','var(\omega_y)','var(\omega_z)')
grid on
% Set y-axis limits
ylim([0, 5*mean([mean(var_wx), mean(var_wy), mean(var_wz)])]);

%% Extract states
phi   = states(:,1);
theta = states(:,2);
psi   = states(:,3);

phi_meas   = meas(:,1);
theta_meas = meas(:,2);
psi_meas   = meas(:,3);

time = Dataset.Time_s_;

%% Plot estimated Euler angles vs measurements (GT)
figure;
subplot(3,1,1)
plot(time, phi, 'r', 'LineWidth', 2); hold on;
plot(time, phi_meas, 'm--', 'LineWidth', 1.5); % GT in magenta
xlabel('Time [s]')
ylabel('\phi [rad]')
title('Roll (\phi) - Estimated vs Measured')
legend('Estimated','Measured (GT)')
grid on

subplot(3,1,2)
plot(time, theta, 'g', 'LineWidth', 2); hold on;
plot(time, theta_meas, 'c--', 'LineWidth', 1.5); % GT in cyan
xlabel('Time [s]')
ylabel('\theta [rad]')
title('Pitch (\theta) - Estimated vs Measured')
legend('Estimated','Measured (GT)')
grid on

subplot(3,1,3)
plot(time, psi, 'b', 'LineWidth', 2); hold on;
plot(time, psi_meas, 'k--', 'LineWidth', 1.5); % GT in black
xlabel('Time [s]')
ylabel('\psi [rad]')
title('Yaw (\psi) - Estimated vs Measured')
legend('Estimated','Measured (GT)')
grid on
