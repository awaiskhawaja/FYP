clc; clear variables;

matpower = case39bus(); % Load MATPOWER data

% Impedance loads
%{
0 = None
1 = Constant Z (1.1 = Constant Z current)
2 = Constant P
%}
load_type = 1;

%% Fault conditions
fault_bus = 20; % Bus that fault occurs at
Z_f = 0; % Fault impedance

%% Setup devices

% 's' = synchronous generator
% 'c' = converter
% ' ' = no current injected

totalBusses = size(matpower.bus, 1);
totalGen = size(matpower.gen, 1);

busses = ones(1, totalBusses) .* ' '; % Devices connected to bus
genNumber = zeros(1, totalBusses); % Generator number connected to bus index
for i = 1:totalGen
    if matpower.gen(i, 22) == 0
        busses(matpower.gen(i, 1)) = 's';
    elseif matpower.gen(i, 22) == 1 || matpower.gen(i, 22) == 2
        busses(matpower.gen(i, 1)) = 'c';
    end

    genNumber(matpower.gen(i, 1)) = i;
end

%% SM

V_s = zeros(1, totalBusses); % Internal voltages for synchronous generators
for i = 1:totalGen
    if matpower.gen(i, 22) == 0
        V_s(matpower.gen(i, 1)) = matpower.gen(i, 6);
    end
end

x_d = zeros(1, totalBusses); % Transient reactance for synchronous generators
for i = 1:totalGen
    if matpower.gen(i, 22) == 0
        
        x_d(matpower.gen(i, 1)) = ((0.2 * 100) / (matpower.gen(i, 9) / 0.8)) * 1j;

    end
end

%% IBR

k = ones(1, totalBusses) .* 2; % Voltage drop control

% Power references for VSCs
P_ref = zeros(1, totalBusses);
Q_ref = zeros(1, totalBusses);
for i = 1:totalGen
    if matpower.gen(i, 22) == 1 || matpower.gen(i, 22) == 2

        P_ref(matpower.gen(i, 1)) = matpower.pfResults.gen(i, 2) / matpower.gen(i, 7);
        Q_ref(matpower.gen(i, 1)) = matpower.pfResults.gen(i, 3) / matpower.gen(i, 7);

    end
end

I_rated = zeros(1, totalBusses);

for i = 1:totalGen

    % S / V
    %S = sqrt(matpower.gen(i, 9)^2 + matpower.gen(i, 4)^2); % 534.7457 MVA
    S = matpower.gen(i, 9) / 0.8;
    b = matpower.bus(matpower.gen(i, 1), :);
    I_rated(matpower.gen(i, 1)) = S/100;

end

overCurrent = 1.5;
I_max = I_rated .* overCurrent;

%% Setup connections

% Generate Y bus
Y_bus = full(matpower.y);

z_org = inv(Y_bus);

% Add subtransient reactances
for i = 1:totalBusses
    if busses(i) == 's' || busses(i) == 'b'
        Y_bus(i, i) = Y_bus(i, i) + (1 / x_d(i));
    end
end

if load_type == 1
    for i = 1:totalBusses
        if matpower.bus(i, 3) ~= 0 || matpower.bus(i, 4) ~= 0
            p_load = (matpower.bus(i, 3) + 1j*matpower.bus(i, 4)) / matpower.baseMVA;
            [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(i, 9)), matpower.pfResults.bus(i, 8));
            v_load = (x + 1j*y);
            i_load = conj(p_load / v_load);
            z_load = v_load / i_load;
            Y_bus(i, i) = Y_bus(i, i) + (1 / z_load);
        end
    end
end

Z_bus = inv(Y_bus);

I_inj_s = zeros(1, totalBusses)'; % SM
I_inj_c = zeros(1, totalBusses)'; % Converters
I_inj_l = zeros(1, totalBusses)'; % Loads

for i = 1:totalBusses

    if busses(i) == 's' || busses(i) == 'b'
        I_inj_s(i) = V_s(i) / x_d(i); % Current injected by s
    end

    if busses(i) == 'c' || busses(i) == 'b'
        % Current injected by c under steady-state pre-fault conditions
        I_inj_c(i) = P_ref(i) + 1j * Q_ref(i);
    end

end

%% Fault calculation
V = zeros(totalBusses, 1);
for i = 1:totalBusses
    [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(i, 9)), matpower.pfResults.bus(i, 8));
    V(i) = (x + 1j*y);
end

I_f_SM = V(fault_bus) / Z_bus(fault_bus, fault_bus);

delta_V = Z_bus(:, fault_bus) * I_f_SM;
V_f = V - delta_V;
I_f_IBR = 0;

for i = 1:totalBusses

    if busses(i) == 'c'
        
        pre_S = P_ref(i) + 1j*Q_ref(i);
        pre_I = conj(pre_S / V(i));
        
        dv = abs(V_f(i)) - abs(V(i));
        delta_I = k(i) * dv * I_rated(i) * 1j;
        new_I = pre_I + delta_I;

        if abs(new_I) > I_max(i) % Exceeding overcurrent capacity

            if abs(imag(new_I)) > I_max(i) % Full Q
                new_I = 1j * I_max(i) * sign(imag(new_I));
                disp("full")
            else % Partial P
                I_Q = imag(new_I); % Retain same Q
                I_P = sign(real(new_I)) * sqrt(I_max(i)^2 - I_Q^2); % Lower P to meet limit
                new_I = I_P + 1j*I_Q;
            end

        end
        
        delta_I = new_I - pre_I;
        new_I_f = (delta_I * Z_bus(i, fault_bus)) / Z_bus(fault_bus, fault_bus);
        I_f_IBR = I_f_IBR + new_I_f;

    end

end

if load_type == 2
    
    p0_load = (matpower.bus(:, 3) + 1j.*matpower.bus(:, 4)) ./ matpower.baseMVA;
    [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(:, 9)), matpower.pfResults.bus(:, 8));
    v0_load = (x + 1j.*y);
    i0_load = -conj(p0_load ./ v0_load);

    i_f_load = -conj(p0_load ./ V_f);
    i_f_load(fault_bus) = 0;

    delta_I = i_f_load - i0_load;
    new_I_f = delta_I .* Z_bus(:, fault_bus);
    new_I_f = sum(new_I_f) / Z_bus(fault_bus, fault_bus);
    I_f_IBR = I_f_IBR + new_I_f;

end

I_f = I_f_SM + I_f_IBR;

endTime = posixtime(datetime('now'));

for i = 1:2:totalBusses
    if i ~= totalBusses
        disp(['V_', num2str(i), ': ', num2str(abs(V_f(i)), 3), ' ∠ ', num2str(rad2deg(angle(V_f(i))), 3), '°', char(9), char(9), char(9),'V_', num2str(i+1), ': ', num2str(abs(V_f(i+1)), 3), ' ∠ ', num2str(rad2deg(angle(V_f(i+1))), 3), '°']);
    else
        disp(['V_', num2str(i), ': ', num2str(abs(V_f(i)), 3), ' ∠ ', num2str(rad2deg(angle(V_f(i))), 3), '°']);
    end
end

disp(' ');

disp(['I_f: ', num2str(abs(I_f), 3), ' ∠ ', num2str(rad2deg(angle(I_f)), 3), '°']);