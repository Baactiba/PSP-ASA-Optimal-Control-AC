function [u_opt, x_opt] = trajOptimizer3(x_initial)
    import casadi.*
    %input current state, output control vector and potential state
    % Define constants
    g = 9.8; % Gravity in m/s^2
    m = 100000; % Mass in kg
    min_thrust = 880 * 1000;  % N
    max_thrust = 2210 * 1000; % N 
    
    length_rod = 60; % Length in meters
    width = 10; % Width 
    
    % Moment of inertia for a uniform density rod
    I = (1/12) * m * length_rod^2; 
    
    % Gimbal angle limits
    max_gimbal = deg2rad(20);
    min_gimbal = -max_gimbal;
    
    % Define the optimization problem
    opti = Opti();
    
    % Set the number of steps and the timestep (dt)
    steps = 400;
    t_step = 0.04;
    
    % Generate the array of state and control vectors
    
    % States: x, x_dot, y, y_dot, theta, theta_dot
    x = opti.variable(steps, 6);  % steps x 6 matrix
    % Controls: thrust (percent), thrust_angle (rad)
    u = opti.variable(steps, 2);
    
    % Initial and final conditions
    % 0, 0, 1000, -80, -pi/2, 0
    
    opti.subject_to(x(1, :) == [x_initial(1), x_initial(2), x_initial(3), x_initial(4), x_initial(5), x_initial(6)]); % Initial state
    opti.subject_to(x(steps, :) == [0, 0, 0, 0, 0, 0]); % Final state
    
    % Cost function to minimize effort and angular velocity
    
    cost = sum(u(:,1).^2) + sum(u(:,2).^2) + 2 * sum(x(:,6).^2);
    opti.minimize(cost);
    
    %constraints
    for i = 1:(steps-1)
        % Current state
        x_current = x(i, :)';
        u_current = u(i, :)';
        
        % Extract state variables
        pos_x = x_current(1);
        vel_x = x_current(2);
        pos_y = x_current(3);
        vel_y = x_current(4);
        theta = x_current(5);
        omega = x_current(6);
        
        % Extract control variables
        thrust_percent = u_current(1);
        thrust_angle = u_current(2);
        
        % forces
        F_x = max_thrust * thrust_percent * sin(thrust_angle + theta);
        F_y = max_thrust * thrust_percent * cos(thrust_angle + theta);
        
        % torque
        T = - (length_rod / 2) * max_thrust * thrust_percent * sin(thrust_angle);
        
        % accelerations
        acc_x = F_x / m;
        acc_y = (F_y / m) - g;
        alpha = T / I;  %angular
        
        % Define the state derivatives
        x_dot = [vel_x; acc_x; vel_y; acc_y; omega; alpha];
        
        % Euler integration for dynamics constraints
        x_next = x_current + x_dot * t_step;
        
        % Impose the dynamics constraint
        opti.subject_to(x(i+1, :)' == x_next);
    end
    
    % Thrust percentage bounds
    opti.subject_to(u(:,1) >= 0.4);
    opti.subject_to(u(:,1) <= 1);
    
    % Thrust angle bounds
    opti.subject_to(u(:,2) >= min_gimbal);
    opti.subject_to(u(:,2) <= max_gimbal);
    
    
    % Solver options
    opti.solver('ipopt');
    
    % Initial guess 
    opti.set_initial(x, 0);
    opti.set_initial(u, 0);
    
    % Solve the optimization problem
    try
        sol = opti.solve();
    catch ME
        disp('An error occurred during optimization:');
        disp(ME.message);
        return;
    end
    
    
    u_opt = sol.value(u);
    x_opt = sol.value(x);
    
    %Plots
    
    figure('Name', 'Optimization Results', 'NumberTitle', 'off', 'Color', 'w');
    
    % Plot state variables
    subplot(2,1,1);
    hold on;
    plot(x_opt(:,1), 'LineWidth', 1.5, 'DisplayName', 'x (Position)');
    plot(x_opt(:,2), 'LineWidth', 1.5, 'DisplayName', 'x\_dot (Velocity)');
    plot(x_opt(:,3), 'LineWidth', 1.5, 'DisplayName', 'y (Position)');
    plot(x_opt(:,4), 'LineWidth', 1.5, 'DisplayName', 'y\_dot (Velocity)');
    plot(x_opt(:,5), 'LineWidth', 1.5, 'DisplayName', 'theta (Angle)');
    plot(x_opt(:,6), 'LineWidth', 1.5, 'DisplayName', 'theta\_dot (Angular Velocity)');
    hold off;
    legend('Location', 'best');
    xlabel('Time Step');
    ylabel('State Values');
    title('State Variables');
    grid on;
    
    % control inputs
    subplot(2,1,2); 
    hold on;
    plot(u_opt(:,1), 'LineWidth', 1.5, 'DisplayName', 'Thrust %');
    plot(u_opt(:,2), 'LineWidth', 1.5, 'DisplayName', 'Thrust Angle (rad)');
    hold off;
    legend('Location', 'best');
    xlabel('Time Step');
    ylabel('Control Inputs');
    title('Control Inputs');
    grid on;
    
    % times
    
    final_time_step = t_step;
    duration = t_step * steps;
    
    fprintf('Final Time Step: %.4f seconds\n', final_time_step);
    fprintf('Total Duration: %.4f seconds\n', duration);

    fprintf("Complete state matrix: \n______________________________\n")
    disp(x_opt)
    fprintf("Complete control matrix: \n______________________________\n")
end