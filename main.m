clc
clear all
global par

par.alpha = 15; % Linear Class-K coefficient
par.alpha_b = 1; % Linear Class-K coefficient for Backup constraint
par.T = 8; % Prediction time
par.Nc = 100; % Number of constraints
par.N = 6; % Number of mixed state-input constraints
par.theta_i = linspace(0,par.T,par.Nc); % Time horizon discretization
par.n = 2; % Dimension of vector x

par.phi_max = 1.75; % Maximal angle
par.u_max = 1.2; % Maximal control input
par.u_min = -1.1; % Minimal control input
par.P_max = 0.2; % Maximal power
par.P_min = -0.7; % Minimal power

par.c = par.phi_max^2; % Backup set size parameter
par.mu = 0.001; % Backup set parameter
par.K = 0.7; % Backup control gain
par.K_b = 0.5*(par.K-sqrt(par.K^2 - 8*par.mu));
par.A_cl = [0, 1; 0, -par.K];

options = odeset('AbsTol',1e-7,'RelTol',1e-8);
initial_cond = [-1.7, 1.9];
[t_b,x_b] = ode45(@main_sim_b, [0 10], initial_cond, options); % Backup simulation
[t_h,x_h] = ode45(@main_sim_h, [0 10], initial_cond, options); % HOCBF simulation

% Backup post-process
input_b = zeros(length(t_b), 1);
H_values_b = zeros(length(t_b), par.N);
for k = 1:length(t_b)
    curr_x_b = x_b(k,(1:2))';
    input_b(k) = k_safe(curr_x_b);   
    curr_u_b = input_b(k);    
    value_a_b = a(curr_x_b);
    value_b_b = b(curr_x_b);
    H_values_b(k,:) = value_b_b - (value_a_b.*curr_u_b);
end

% HOCBF post-process
input_h = zeros(length(t_h), 1);
H_values_h = zeros(length(t_h), par.N);
for k = 1:length(t_h)
    curr_x_h = x_h(k,(1:2))';
    input_h(k) = k_safe(curr_x_h);   
    curr_u_h = input_h(k);    
    value_a_h = a(curr_x_h);
    value_b_h = b(curr_x_h);
    H_values_h(k,:) = value_b_h - (value_a_h.*curr_u_h);
end
%% Plots
figure(1)
Phi_stack = zeros(2, 2, length(par.theta_i));
for j = 1:length(par.theta_i)
    Phi_stack(:, :, j) = expm(par.A_cl * par.theta_i(j));
end
Invariant_set = @(X, Y) arrayfun(@(x, y) min_calculation(x, y, Phi_stack, @k_b, @a, @b, @h_b), X, Y);
KB_handle = @(X, Y) arrayfun(@(x, y) k_b([x, y]), X, Y);
subplot(3,2,[1 3 5])
fimplicit(Invariant_set, 'Color', 'g', 'LineWidth', 1.2);
hold on
plot(x_b(:,1), x_b(:,2),'LineWidth',2, 'Color','b')
hold on
plot(x_h(:,1), x_h(:,2),'LineWidth',2, 'Color', [1 0.5 0], 'LineStyle','--')
hold on
xlim([-1.8 1.8])
ylim([-2 2])
xlabel('angle, $x_1$ (rad)','interpreter','latex')
ylabel('ang. velocity, $x_2$ (rad/s)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
set(gca,'fontsize', 14) 
hold on
fimplicit(@(x,y) min(-2*x.*y + par.K*(par.phi_max^2 - x.^2),(par.phi_max^2 - x.^2)), 'Color','k') % HOCBF

subplot(3,2,2)
plot(t_b,x_b(:,1), 'LineWidth',2, 'Color','b')
hold on
plot(t_h,x_h(:,1), 'LineWidth',2, 'Color', [1 0.5 0], 'LineStyle','--')
hold on
yline(1.75,'LineStyle','--')
yline(-1.75,'LineStyle', '--')
ylim(1.1*1.75 * [-1 1])
xlim([0 10])
ylabel('angle, $x_1$ (rad)','interpreter','latex')
xlabel('time, $t$ (s)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
set(gca,'fontsize', 14) 

subplot(3,2,4)
plot(t_b,input_b, 'LineWidth',2, 'Color','b')
hold on
plot(t_h,input_h, 'LineWidth',2, 'Color', [1 0.5 0], 'LineStyle','--')
hold on
yline(par.u_max,'LineStyle','--')
yline(par.u_min,'LineStyle','--')
ylim([-1.5 1.5])
xlim([0 10])
xlabel('time, $t$ (s)','interpreter','latex')
ylabel('input, $u$ (Nm)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
set(gca,'fontsize', 14) 
hold on

subplot(3,2,6)
plot(t_b,x_b(:,2).*input_b, 'LineWidth',2, 'Color','b')
hold on
plot(t_h,x_h(:,2).*input_h, 'LineWidth',2, 'Color', [1 0.5 0], 'LineStyle','--')
hold on
yline(par.P_max,'LineStyle','--')
yline(par.P_min,'LineStyle','--')
ylim([par.P_min-0.1 par.P_max+0.1])
xlim([0 10])
xlabel('time, $t$ (s)','interpreter','latex')
ylabel('power, $P$ (W)','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
set(gca,'fontsize', 14) 
hold on
%% Functions
function xdot = main_sim_b(t,x)
disp(t)
xdot = f(x) + g(x)*k_safe(x);
end

function xdot = main_sim_h(t,x)
disp(t)
xdot = f(x) + g(x)*k_h(x);
end

function outCBF = k_h(x)
global par
h = -2*x(1)*x(2)+par.K*(par.phi_max^2 - x(1)^2);
gradh = [-2*x(2)-2*par.K*x(1), -2*x(1)];
H = eye(1);
F = [];
Aeq = [];
beq = [];
x0 = 0;
A = [];
B = [];
aa = a(x);
bb = b(x);
A = [A; -gradh*g(x); aa(3:6)];
B = [B; par.alpha*h+gradh*f(x)+gradh*g(x)*k_d(x); bb(3:6)];
options = optimoptions('quadprog','Algorithm','active-set');
outCBF = quadprog(H,F,A,B,Aeq,beq,[],[],x0,options) + k_d(x);
end

function outQP = k_safe(x)
global par
H = eye(1);
F = [];
Aeq = [];
beq = [];
x0 = 0;
A = [];
B = [];
for ii = 1:par.Nc
    phi_i = expm(par.A_cl*par.theta_i(ii))*x;
    Q_i = expm(par.A_cl*par.theta_i(ii));
    [h_phi_i, gradh_phi_i] = CBF(phi_i);
    for j = 1:par.N
        A = [A; -gradh_phi_i(j,:)*Q_i*g(x);];  
        B = [B; par.alpha*h_phi_i(j)+gradh_phi_i(j,:)*Q_i*f(x)+gradh_phi_i(j,:)*Q_i*g(x)*k_d(x)];
    end
end
phi_T = expm(par.A_cl*par.theta_i(end))*x;
Q_T = expm(par.A_cl*par.theta_i(end));
[h_b_phi_T,gradhb_phi_T] = h_b(phi_T);
A = [A; -gradhb_phi_T*Q_T*g(x); a(x)];
B = [B; par.alpha_b*h_b_phi_T+gradhb_phi_T*Q_T*f(x)+gradhb_phi_T*Q_T*g(x)*k_d(x); b(x)];
options = optimoptions('quadprog','Algorithm','active-set');
outQP = quadprog(H,F,A,B,Aeq,beq,[],[],x0,options) + k_d(x);
end

function min_margin = min_calculation(x_val, y_val, Phi_stack, k_b_fun, a_fun, b_fun, h_b_fun)
    z0 = [x_val; y_val];
    num_steps = size(Phi_stack, 3);
    % Initialize large value
    min_margin = inf; 
    for k = 1:num_steps
        % 1. Evolve state: x(t) = Phi(t) * x(0)
        phi_t = Phi_stack(:, :, k) * z0;
        % 2. Evaluate running constraint
        u_b = k_b_fun(phi_t'); 
        a_vec = a_fun(phi_t');
        b_vec = b_fun(phi_t');
        % 3. Compute b - a*u >= 0
        margins = b_vec - (a_vec .* u_b);
        current_worst = min(margins);
        if current_worst < min_margin
            min_margin = current_worst;
        end
    end
    % Terminal constraint at time T
    phi_T = Phi_stack(:, :, end) * z0;
    terminal_margin = h_b_fun(phi_T');
    if terminal_margin < min_margin
        min_margin = terminal_margin;
    end
end