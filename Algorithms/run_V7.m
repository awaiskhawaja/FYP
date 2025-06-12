function results = run_V7(matpower, fault_bus, Z_f, load_type, FRT, overCurrent)
    
    %% Loads
    %{
    0 = None
    1 = Constant Z (1.1 = Constant Z current)
    2 = Constant P
    %}
    %load_type = 1;

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
            
            x_d(matpower.gen(i, 1)) = ((0.2 * 100) / (matpower.gen(i, 9) / 0.8))  * 1j;

        end
    end

    %% IBR

    k = ones(1, totalBusses) .* FRT; % Voltage drop control
    
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
    
    %overCurrent = 1.1;
    I_max = I_rated .* overCurrent;

    %% Setup connections

    % Generate Y bus
    Y_bus = full(matpower.y);
    
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

    %% Pre-Fault Setup

    V = [];
    I_f = 0;
    delta_V = [];
    saturation_state = zeros(1, totalBusses);
    automatic_state = 1;
    
    tolerance = 0.001;
    
    running = 1;
    m = 0; % Iteration number
    m_max = 10000;
    
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
    I_IBR_pf = I_inj_c;

    % SM contribution to fault
    [x, y] = pol2cart(deg2rad(matpower.pfResults.bus(fault_bus, 9)), matpower.pfResults.bus(fault_bus, 8));
    v_fault_bus = (x + 1j*y);
    I_RP = v_fault_bus / (Z_bus(fault_bus, fault_bus) + Z_f);
    I_f = I_RP;
    
    V = Z_bus * (I_inj_s + I_inj_c + I_inj_l); % Voltage caused by injected current
    delta_V = Z_bus(:, fault_bus) * I_RP; % Voltage drop caused by fault
    V = V - delta_V; % New voltage at busses
    old_V = V;

    %% Iteration

    while running
        
        % Optional, increase tolerance if non-convergence
        % if m > 200
        %     if mod(m, 200) == 0
        %         tolerance = tolerance * 2;
        %     end
        % end

        % Update converter current references
        all_converged = 1;
        for i = 1:totalBusses
    
            if busses(i) == 'c' || busses(i) == 'b'
                
                % After 5000 iterations, force it to end
                if m > 5000
                    V = old_V;
                end
                
                 % Update load currents
                if load_type == 2
                    I_inj_l = find_I_load(matpower, V, fault_bus);
                end
    
                if load_type == 1.1
                    I_inj_l = find_I_load_Z(matpower, V, fault_bus);
                end

                V = Z_bus * (I_inj_s + I_inj_c + I_inj_l); % Voltage caused by injected current

                if m > 5000
                    V = old_V;
                end

                % Fault current
                I_CDP = 0; % IBR
                for b = 1:totalBusses
                    if busses(b) == 'c'
                        delta_I = I_inj_c(b) - I_IBR_pf(b);
                        I_CDP = I_CDP + (delta_I * Z_bus(fault_bus, b));
                    end
                end
                I_CDP = I_CDP / (Z_bus(fault_bus, fault_bus) + Z_f);

                I_L = 0; % Load
                for b = 1:totalBusses
                    if abs(I_load_pf(b)) > 0
                        delta_I = I_inj_l(b) - I_load_pf(b);
                        I_L = I_L + (delta_I * Z_bus(fault_bus, b));
                    end
                end
                I_L = I_L / (Z_bus(fault_bus, fault_bus) + Z_f);
                I_f = I_RP + I_CDP + I_L;

                delta_V = Z_bus(:, fault_bus) * I_f; % Voltage drop caused by fault
                V = V - delta_V; % New voltage at busses
                
                % After 2000 iterations, limit rate of change of voltage
                % between iterations to help convergence
                if m > 2000
                    a = 0.9 - (m/1000);
                    %a = 0.5;
                    if a < 0.05
                        a = 0.05;
                    end
                    V = (V * a) + (old_V * (1 - a));
                end

                if m > 5000
                    V = old_V;
                end

                % Update Q reference based on voltage drop
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

        if m > 500
            a = 0.5 - (m/10000);
            %a = 0.5;
            if a < 0.05
                a = 0.05;
            end
            V = (V * a) + (old_V * (1 - a));
        end
    
        old_V = V;
    
        if all_converged
            running = 0; % Break loop
        end
    
        m = m + 1;
        
        % Debugging
        % store_V(m) = abs(V(30));
        %store_V(m) = abs(I_f);
        %store_V(m) = abs(I_inj_c(37));
    
        if m >= m_max % Iteration limit if not converged
            running = 0;
        end
        
    end

    results.I_f = abs(I_f); % Overall fault current magnitude
    results.m = m; % Total iterations
    results.tolerance = tolerance; % Tolerance
    results.V = V; % Bus voltage (complex)
    results.I_inj_c = I_inj_c; % IBR current injections (complex)
    results.I_inj_c_PU = abs(I_inj_c) ./ I_rated'; % IBR current injection magnitudes per unit

end