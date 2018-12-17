%% Task 4 Problem Statement

%% Constants/Assumptions

mdot_air = 6.0;
M_air = 28.97;
P1 = 101; %[kPa]
T1 = 298; %[K], 25oC
P2 = 500; %[kPa]
T4 = 1600; %[K]
g = 0.25;
Qs = 0;
eff_comp = 0.85;
eff_turb = 0.85;
e_regen = 0.75;

h1 = 0; %[kJ/kmol] - find enthalpies at other states by determining change in enthalpy
R = 8.314; %[kJ/(kmol*K)]

a_stoich = 4.76*(2+3*g);
a = a_stoich;

for z = 1:2
    %% (i) Compressor
    [h2, T2, W_req] = Task2(P2, eff_comp);

    %% (ii)
    % Implement Task3 where necessary
    
    %% (iii & iv): Turbine Expansion Process
    error = [1000 1000 1000];
    T5_guess = [500 0 2000];
    tol = 0.2;

    i=1;
    while(abs(error(2))>tol)
        mid = (T5_guess(1)+T5_guess(3))/2;
        T5_guess(2) = mid;

        Tavg = 0.5*(T4+T5_guess);
        cp_prod_turbine = Task3(a,g,Tavg);
        cv_prod_turbine = cp_prod_turbine - R;
        k = cp_prod_turbine./cv_prod_turbine;
        T5_s = T4*(P1/P2).^((k-1)./k);
        T5_calc = eff_turb.*(T5_s - T4) + T4;

        error = T5_calc - T5_guess;

        %Tighten the T5_guess bounds
        if error(1)*error(2)<0
            T5_guess(3) = T5_guess(2);
        else
            T5_guess(1) = T5_guess(2);
        end

        T5_result(i) = T5_guess(2);
        T5_constant(i) = T5_calc(2);
        cp_prod_turbine_result(i) = cp_prod_turbine(2);
        i=i+1;
    end
    T5 = T5_result(end);
    cp_prod_turbine_out = cp_prod_turbine_result(end);

    
    %% (v): Regenerator
    
    Tavg_regen = 0.5*(T2 + T5_result(end));
    [cp_air_regen, cp_prod_regen] = Task3(a,g,Tavg_regen);
    
    M_N2 = 28.01;
    M_O2 = 32.00;
    M_H2O = 18.02;
    M_CO2 = 44.01;
    
    n_prod = 3.76/4.76*a + (a/4.76-3*g-2) + (2+2*g) + (1+2*g);
    m_prod = M_N2*3.76/4.76*a + M_O2*(a/4.76-3*g-2) + M_H2O*(2+2*g) + ...
        M_CO2*(1+2*g);
    
    ndot_prod = mdot_air/(m_prod/n_prod); % [kmol/s]
    ndot_air = mdot_air/M_air; % [kmol/s]
    
    n_cp_air = ndot_air*cp_air_regen;
    n_cp_prod = ndot_prod*cp_prod_regen;
    n_cp_min = min([n_cp_air n_cp_prod]);

    syms T2_r
    eqn1 = e_regen == (n_cp_air)*(T2_r - T2)/(n_cp_min*(T5-T2));
    T2_r_ans = double(solve(eqn1)); % [K]

    %% (vi): Find T3
    
    [f, cp_air] = Task3(a,g,T2_r_ans);
    T3 = Qs/(ndot_air*cp_air) + T2_r_ans;

    % can only evaluate cp_air at T2_r_ans if T2_r_ans = T3; 
    % this is the case for Task 4; however, for Task 5, T2_r_ans != T3;
    % For Task 5, cp_air must be evaluated at Tavg of T2_r and T3.
    
    %% (vii) Compute a
    % Eq. 3
    a_corrected = Task4supp(g, T3, T4);
    a_result(z) = a_corrected;
end
a_result = a_result(1);

%% Performance Calculations
W_out = ndot_prod*cp_prod_turbine_out*(T4-T5);
%Tavg_burner = 0.5*(T4-T3);
Tavg_burner = T4;
cp_prod_burner = Task3(a_result,g,Tavg_burner);

Qburner = ndot_prod*cp_prod_burner*(T4-T3); %Q for combustor
power_out = W_out-W_req; %[W]
eff_sys = power_out/(Qs + Qburner);

eff_class = (eff_turb*(1-(P1/P2)^((k-1)/k)) - (1/eff_comp*T1/T4*((P2/P1)^((k-1)/k)-1)))/...
    (1-e_regen*(1-eff_turb+eff_turb*(P1/P2)^((k-1)/k)) - ...
     (1-e_regen)*(T1/T4)*(1-(1/eff_comp) + (1/eff_comp)*(P2/P1)^((k-1)/k)));

%% Plot
% plot(T5_result);
% hold on;
% plot(T5_constant);
% xlabel("# of Iterations");
% ylabel("Temperature [K]");
% legend("T5 Guess","T5 Calculated","Location","east");