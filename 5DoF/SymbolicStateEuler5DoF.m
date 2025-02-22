% constants (same as 3DoF)
g = 9.81; % [m / s2]

syms mass L I; 
c = [mass; L; I];

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
theta = sym("theta", [3;1]);
w = sym("w", [2;1]);

x = [r; v; theta; w];

% controls

T = sym("T", [3;1]);;

u = [norm_T];

% calculate theta dot
b_inverse =  (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3))
     0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3))
     sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];
zero_w = [0; w];
thetadot = b_inverse * zero_w;

% r dot
r_dot = v;

Rxtheta1 = [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
Rxtheta2 = [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
Rxtheta3 = [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
C_be = Rxtheta3 * Rxtheta2 * Rxtheta1;
C_eb = C_be.';
T_e = C_eb * T;
% v dot
v_dot = T_e / mass - [0; 0; g];

% w dot
M = cross([-L; 0; 0], T);
w_dot = M / I;

x_dot = [r_dot; v_dot; thetadot; w_dot(2:3)];
j_a = jacobian(x_dot, x);
j_b = jacobian(x_dot, u);

vars = [x; u; mass; L; I];


matlabFunction(x_dot, "File", "5DoF/SymDynamics5DoF", "Vars", {[x], [u], [c]})
matlabFunction(j_a, "File", "5DoF/SymDynamics5DoF_j_a", "Vars", {[x], [u], [c]})
matlabFunction(j_b, "File", "5DoF/SymDynamics5DoF_j_b", "Vars", {[x], [u], [c]})


% Create equations of motion block for Simulink model
open("5DoF/EoM_Euler5DoF.slx")
matlabFunctionBlock('SymDynamicsEuler5DoFBlock',x_dot,'Vars',{[x], [u], [c]});
close_system("5DoF/EoM_Euler5DoF.slx", 1)

%open("JacobianX_5DoF.slx")
%matlabFunctionBlock("JacobianX_5DoF/SymXJacobianEuler5DoF",j_a,"Vars",vars);
%close_system("JacobianX_5DoF.slx", 1)