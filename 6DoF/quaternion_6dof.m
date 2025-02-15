% constants (same as 3DoF)
syms mass L I1 I2 I3 g; 
c = [g; mass; L; I1; I2; I3]

% states

r = sym("r", [3;1])
v = sym("v", [3;1])
q = sym("q", [4;1])
w = sym("w", [3;1])

% controls

x = [r; v; q; w]

T = sym("T", [3;1])
gamma = sym("gamma", 1)

u = [T; gamma]

% calculate q dot
matrix_with_w =  [0 -w(1) -w(2) -w(3); w(1) 0 w(3) -w(2);
    w(2) -w(3) 0 w(1); w(3) w(2) -w(1) 0]

matrix_with_q =  [q(1) -q(2) -q(3) -q(4); q(2) q(1) q(4) -q(3);
    q(3) -q(4) q(1) q(2); q(4) q(3) -q(2) q(1)]

q_dot_fixed = 0.5 * (matrix_with_w * q)
q_dot_body = 0.5 * (matrix_with_q * [0;  w]) % use this one for now

% x dot

r_dot = v

% v dot
rotation_matrix = [q(1)*q(1)+q(2)*q(2)-q(3)*q(3)-q(4)*q(4)   2*(q(2)*q(3) + q(4)*q(1))  2*(q(2)*q(4) - q(3)*q(1));
                   2*(q(2)*q(3) - q(4)*q(1))  q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)   2*(q(3)*q(4) + q(2)*q(1));
                   2*(q(2)*q(4) + q(3)*q(1))  2*(q(3)*q(4) - q(2)*q(1))  q(1)*q(1) + q(4)*q(4) - q(3)*q(3) - q(2)*q(2)]
C_eb = rotation_matrix.'
v_dot = rotation_matrix * T / c(2) * cos(u(4)) + [0; 0; -g]

% w dot

mag = sqrt(sum(T.^2))
M = sin(u(4))*sqrt(sum(T.^2))*[1; 0; 0] - cross([c(3); 0; 0], cos(u(4)) * T)
w_dot = (M + (c([5; 6; 4]) - c([6; 4; 5])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ c([4; 5; 6]) 



x_dot = [r_dot; v_dot; q_dot_body; w_dot]
j_a = jacobian(x_dot, x)
j_b = jacobian(x_dot, u)


xdot_func6q = matlabFunction(x_dot);
j_a_func6q = matlabFunction(j_a);
j_b_func6q = matlabFunction(j_b);

save("quaternion_6dof.m", "xdot_func6q", "j_a_func6q", "j_b_func6q")
save("SymDynamics6DoFQ.m", "xdot_func6q")
save("SymUJacobian6DoFQ.m", "j_b_func6q")
save("SymXJacobian6DoFQ.m", "j_a_func6q")