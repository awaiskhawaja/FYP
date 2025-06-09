clc; clear variables;

matpower = case39bus_3(1, 1); % Load MATPOWER data

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
% 'b' = both
% ' ' = no current injected

totalBusses = size(matpower.bus, 1);
totalGen = size(matpower.gen, 1);

busses = ones(1, totalBusses) .* ' '; % Devices connected to bus
genNumber = zeros(1, totalBusses);
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

        %I_sm = (matpower.gen(i, 2) + 1j * matpower.gen(i, 3)) / 100;
        %x_d(matpower.gen(i, 1)) = (matpower.gen(i, 6) - 1) / I_sm;
    end
end

%% IBR

k = ones(1, totalBusses) .* 2; % Voltage drop control

% Power references for VSCs
P_ref = zeros(1, totalBusses);
Q_ref = zeros(1, totalBusses);
for i = 1:totalGen
    if matpower.gen(i, 22) == 1 || matpower.gen(i, 22) == 2

        %P_ref(matpower.gen(i, 1)) = matpower.gen(i, 2) / matpower.gen(i, 7);
        %Q_ref(matpower.gen(i, 1)) = matpower.gen(i, 3) / matpower.gen(i, 7);

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
Y_bus = full(matpower.y); % Only need first 10x10 values (since bus 11 is a "virtual bus")

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

%% Initial calculations

I_inj_s = zeros(1, totalBusses)'; % SM
I_inj_c = zeros(1, totalBusses)'; % Converters
I_inj_l = zeros(1, totalBusses)'; % Loads

for i = 1:totalBusses

    if busses(i) == 's' || busses(i) == 'b'
        I_inj_s(i) = V_s(i) / x_d(i); % Current injected by s
    end

    if busses(i) == 'c' || busses(i) == 'b'
        % Current injected by c under steady-state pre-fault conditions
        [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(i, 9)), matpower.pfResults.bus(i, 8));
        v0 = (x + 1j*y);
        I_inj_c(i) = conj((P_ref(i) + 1j * Q_ref(i)) / v0);
    end

end

%% Iteration

V = [];
I_f = 0;
I_f_G74 = 0
delta_V = [];
saturation_state = zeros(1, totalBusses);
automatic_state = 1;

tolerance = 0.001;

running = 1;
m = 0; % Iteration number
m_max = 10000;

% Data logging
store_V = zeros(1, m_max);
P_store = zeros(1, totalBusses);
Q_store = zeros(1, totalBusses);
SM_store = zeros(1, totalBusses);

function I_load = find_I_load(matpowerCase, voltages, f_bus)
    p0_load = (matpowerCase.bus(:, 3) + 1j.*matpowerCase.bus(:, 4)) ./ matpowerCase.baseMVA; % Constant power
    I_load = -conj(p0_load ./ voltages); % Current drawn from bus
    I_load(f_bus) = 0;
end

function I_load_Z = find_I_load_Z(matpowerCase, voltages, f_bus)
    p0_load = (matpowerCase.bus(:, 3) + 1j.*matpowerCase.bus(:, 4)) ./ matpowerCase.baseMVA;
    [x, y] = pol2cart(deg2rad(matpowerCase.pfResults.bus(:, 9)), matpowerCase.pfResults.bus(:, 8));
    v0_load = (x + 1j.*y);
    p_new_load = ((abs(voltages).^2) ./ (abs(v0_load).^2)) .* p0_load;
    I_load_Z = -conj(p_new_load ./ voltages);
    I_load_Z(f_bus) = 0;
end

% Pre-fault conditions
if load_type == 2
    [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(:, 9)), matpower.pfResults.bus(:, 8));
    v0_load = (x + 1j.*y);
    I_inj_l = find_I_load(matpower, v0_load, fault_bus);
end

if load_type == 1.1
    [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(:, 9)), matpower.pfResults.bus(:, 8));
    v0_load = (x + 1j.*y);
    I_inj_l = find_I_load_Z(matpower, v0_load, fault_bus);
end

I_load_pf = I_inj_l;

[x, y] = pol2cart(deg2rad(matpower.pfResults.bus(:, 9)), matpower.pfResults.bus(:, 8));
v0 = (x + 1j.*y);
I_IBR_pf = conj((P_ref + 1j .* Q_ref).' ./ v0);

% SM contribution to fault
[x, y] = pol2cart(deg2rad(matpower.pfResults.bus(fault_bus, 9)), matpower.pfResults.bus(fault_bus, 8));
v_fault_bus = (x + 1j*y);
I_RP = v_fault_bus / (Z_bus(fault_bus, fault_bus) + Z_f);

V = Z_bus * (I_inj_s + I_inj_c + I_inj_l); % Voltage caused by injected current
I_f = V(fault_bus) / (Z_bus(fault_bus, fault_bus) + Z_f); % Fault current
delta_V = Z_bus(:, fault_bus) * I_RP; % Voltage drop caused by fault
V = V - delta_V; % New voltage at busses

while running

    % Option to increase tolerance if non-convergance. Need to check
    % tolerance/number of iterations to check if it has converged properly
    if m > 200
        if mod(m, 200) == 0
            tolerance = tolerance * 2;
        end
    end
  
    % Update converter current references
    all_converged = 1;
    for i = 1:totalBusses

        if busses(i) == 'c' || busses(i) == 'b'
            
             % Update load currents
            if load_type == 2
                I_inj_l = find_I_load(matpower, V, fault_bus);
            end

            if load_type == 1.1
                I_inj_l = find_I_load_Z(matpower, V, fault_bus);
            end

            V = Z_bus * (I_inj_s + I_inj_c + I_inj_l); % Voltage caused by injected current
            
            % Fault current
            I_f = V(fault_bus) / (Z_bus(fault_bus, fault_bus) + Z_f);
            I_CDP = 0; % IBR fault current contribution
            for j = 1:totalBusses
                if busses(j) == 'c'
                    delta_I = I_inj_c(j) - I_IBR_pf(j);
                    I_CDP = I_CDP + (delta_I * Z_bus(fault_bus, j));
                end
            end
            I_CDP = I_CDP / (Z_bus(fault_bus, fault_bus) + Z_f);

            I_L = 0; % Load currents (modelled in similar way to IBR)
            for j = 1:totalBusses
                if abs(I_load_pf(j)) > 0
                    delta_I = I_inj_l(j) - I_load_pf(j);
                    I_L = I_L + (delta_I * Z_bus(fault_bus, j));
                end
            end
            I_L = I_L / (Z_bus(fault_bus, fault_bus) + Z_f);
            I_f_G74 = I_RP + I_CDP + I_L;

            delta_V = Z_bus(:, fault_bus) * I_f_G74; % Voltage drop caused by fault
            V = V - delta_V; % New voltage at busses
            
            % Update Q reference based on voltage drop
            %Q_temp = abs(V(i)) * (Q_ref(i) + (k(i) * (1 - abs(V(i)))));
            %Q_temp = abs(V(i)) * (Q_ref(i) + (k(i) * I_rated(i) * (1 - abs(V(i)))));
            Q_temp = abs(V(i)) * (Q_ref(i) + (k(i) * I_rated(i) * (abs(matpower.pfResults.bus(i, 8)) - abs(V(i)))));
            P_temp = P_ref(i);
            
            %{
            Saturation states:
            0 = Unsaturated (USS) -> Full P and Q supplied
            1 = Partially Saturated (PSS) -> Full Q, P limited by I_max
            2 = Fully Saturated (FSS) -> P = 0, Q limited by I_max
            %}

            if automatic_state % Determine what saturation state each c is in
                
                % Minimum bus voltage needed for VSC to operate in USS
                USS_lim = sqrt((P_temp^2) + (Q_temp^2)) / I_max(i); 
    
                % Minimum bus voltage needed for VSC to not operate in FSS
                FSS_lim = abs(Q_temp) / I_max(i);
                
                if V(i) >= USS_lim
                    saturation_state(i) = 0;
                elseif V(i) >= FSS_lim
                    saturation_state(i) = 1;
                else
                    saturation_state(i) = 2;
                end

            end
            
            if saturation_state(i) == 1 % PSS
                % Limit P so that Q requirement can still be met, but current doesn't go above I_max
                P_temp = sign(P_temp) * sqrt((abs(V(i)) * I_max(i))^2 - Q_temp^2);
            elseif saturation_state(i) == 2
                P_temp = 0;
            end

            % New current injected by c to meet new power setpoint
            new_i = conj((P_temp + 1j * Q_temp) / V(i));

            if saturation_state(i) == 2 % FSS
                % Limit current to magnitude of I_max
                old_i = new_i;
                new_i = (new_i / abs(new_i)) * I_max(i);
                Q_temp = Q_temp * (abs(new_i) / abs(old_i));
            end
            
            % If difference between iterations is greater than tolerance
            if abs(new_i - I_inj_c(i)) > tolerance
                all_converged = 0; % Set flag to flase
            end

            I_inj_c(i) = new_i;

            P_store(i) = P_temp;
            Q_store(i) = Q_temp;

        end
    end

    if all_converged
        running = 0; % Break loop
    end

    m = m + 1;
    if m >= m_max % Iteration limit if not converged
        running = 0;
    end
    
end


%% Display results

if load_type == 2
    I_inj_l = find_I_load(matpower, V, fault_bus);
end

if load_type == 1.1
    I_inj_l = find_I_load_Z(matpower, V, fault_bus);
end

for i = 1:2:totalBusses
    if i ~= totalBusses
        disp(['V_', num2str(i), ': ', num2str(abs(V(i)), 3), ' ∠ ', num2str(rad2deg(angle(V(i))), 3), '°', char(9), char(9), char(9),'V_', num2str(i+1), ': ', num2str(abs(V(i+1)), 3), ' ∠ ', num2str(rad2deg(angle(V(i+1))), 3), '°']);
    else
        disp(['V_', num2str(i), ': ', num2str(abs(V(i)), 3), ' ∠ ', num2str(rad2deg(angle(V(i))), 3), '°']);
    end
end

disp(' ');

disp(['I_f: ', num2str(abs(I_f_G74), 3), ' ∠ ', num2str(rad2deg(angle(I_f_G74)), 3), '°']);

% Following code is for analysing and displaying other information for
% debugging purposes

% x = 1:totalBusses;
% y = zeros(1, totalBusses);
% sum_delta_I = 0;
% sum_delta_I_G74 = 0;
% sum_delta_I_G74_2 = 0;
% 
% disp(' ');
% for i = 1:totalBusses
% 
%     if busses(i) == 's' || busses(i) == 'c'
% 
%         g = genNumber(i);
% 
%         %pre_S = (matpower.gen(g, 2) + 1j*matpower.gen(g, 3)) / matpower.baseMVA;
%         pre_S = (matpower.pfResults.gen(g, 2) + 1j*matpower.pfResults.gen(g, 3)) / matpower.baseMVA;
%         [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(i, 9)), matpower.pfResults.bus(i, 8));
%         pre_V = (x + 1j*y);
%         pre_I = conj(pre_S / pre_V);
% 
%         if i == 30
%             abs(pre_I)
%         end
% 
%     end
% 
%     if busses(i) == 's'
%         SM_I = (V_s(i) - V(i)) / x_d(i); % SM output current
%         SM_store(i) = V_s(i) * conj(SM_I); % SM output complex power
%         SM_pu = SM_I / I_rated(i); % SM output current in pu
% 
%         delta_I = abs(SM_I) - abs(pre_I);
%         sum_delta_I = sum_delta_I + delta_I;
%         delta_I_pu = delta_I / I_rated(i);
% 
%         disp(['I_', num2str(i), ' (SM): ', num2str(abs(SM_I), 3), ' (', num2str(abs(SM_pu), 3), ' rated)']);
%         disp([9, 9, 'ΔI: ', num2str(delta_I, 3), ' (', num2str(delta_I_pu, 3), ' rated)']);
%         disp(' ');
%         y(i) = abs(SM_I);
% 
%     elseif busses(i) == 'c'
%         IBR_pu = I_inj_c(i) / I_rated(i); % IBR output current in pu
% 
%         delta_I = abs(I_inj_c(i)) - abs(pre_I);
%         sum_delta_I = sum_delta_I + delta_I;
%         delta_I_pu = delta_I / I_rated(i);
% 
%         temp = (I_inj_c(i) - pre_I) * Z_bus(i, fault_bus);
%         %temp = (I_inj_c(i)) * Z_bus(i, fault_bus);
%         sum_delta_I_G74 = sum_delta_I_G74 + temp;
% 
%         % 2nd method
%         dv = abs(matpower.pfResults.bus(i, 8)) - abs(V(i));
%         delta_I = k(i) * dv * I_rated(i) * 1j;
%         new_I = pre_I + delta_I;
% 
%         if abs(new_I) > I_max(i) % Exceeding overcurrent capacity
% 
%             if abs(imag(new_I)) > I_max(i) % Full Q
%                 new_I = 1j * I_max(i) * sign(imag(new_I));
%             else % Partial P
%                 I_Q = imag(new_I); % Retain same Q
%                 I_P = sign(real(new_I)) * sqrt(I_max(i)^2 - I_Q^2); % Lower P to meet limit
%                 new_I = I_P + 1j*I_Q;
%             end
% 
%         end
% 
%         % delta_I = new_I % non-delta
%         delta_I = new_I - pre_I;
%         new_I_f = delta_I * Z_bus(i, fault_bus);
%         sum_delta_I_G74_2 = sum_delta_I_G74_2 + new_I_f;
%         % end 2nd method
% 
%         disp(['I_', num2str(i), ' (IBR): ', num2str(abs(I_inj_c(i)), 3), ' (', num2str(abs(IBR_pu), 3), ' rated)']);
%         disp([9, 9, 'ΔI: ', num2str(delta_I, 3), ' (', num2str(delta_I_pu, 3), ' rated)']);
%         disp(' ');
%         y(i) = abs(I_inj_c(i));
%     end
% 
% end
% 
% % disp(' ');
% % disp(['ΔI (generators): ', num2str(sum_delta_I, 3)]);
% 
% % Change in load currents
% if load_type == 1.1 || load_type == 2
% 
%     p0_load = (matpower.bus(:, 3) + 1j.*matpower.bus(:, 4)) ./ matpower.baseMVA;
%     [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(:, 9)), matpower.pfResults.bus(:, 8));
%     v0_load = (x + 1j.*y);
%     i0_load = -conj(p0_load ./ v0_load);
% 
%     %delta_I = abs(I_inj_l) - abs(i0_load);
%     % sum_delta_load = sum(delta_I);
%     % delta_load_G74 = (I_inj_l - i0_load) .* Z_bus(:,fault_bus);
%     % sum_delta_load_G74 = sum(delta_load_G74);
% 
%     %delta_I = I_inj_l;% non-delta
%     delta_I = I_inj_l - i0_load;
%     sum_delta_load = abs(sum(delta_I));
%     delta_load_G74 = delta_I .* Z_bus(:, fault_bus);
%     sum_delta_load_G74 = sum(delta_load_G74);
% 
%     % disp(['ΔI (loads): ', num2str(sum_delta_load, 3)]);
% else
%     sum_delta_load = 0;
% end
% 
% total_delta_I = sum_delta_I + sum_delta_load;
% % disp(['Total ΔI: ', num2str(total_delta_I, 3)]);
% 
% % disp(' ');
% % disp('G74 Method:');
% 
% [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(fault_bus, 9)), matpower.pfResults.bus(fault_bus, 8));
% v_fault_bus = (x + 1j*y);
% I_RP_after = v_fault_bus / Z_bus(fault_bus, fault_bus);
% I_CDP_after = sum_delta_I_G74 / Z_bus(fault_bus, fault_bus);
% 
% I_CDP_2 = sum_delta_I_G74_2 / Z_bus(fault_bus, fault_bus);
% 
% if load_type == 1.1 || load_type == 2
%     I_L_after = sum(delta_load_G74 ./ Z_bus(fault_bus, fault_bus));
% else
%     I_L_after = 0;
% end
% 
% I_f_G74_after = I_RP_after + I_CDP_after + I_L_after;
% 
% disp(['Fault Current: ', num2str(abs(I_f_G74_after), 3)]);