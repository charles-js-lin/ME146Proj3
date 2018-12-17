%% Task 2 Problem Statement
% Compute exit condition from compressor and work required per kg of air
 function [h2_out, T2_out, W_req] = Task2(P2, eff_comp)
    %% Initialize Constants/Assumptions
%     P2 = 130;
%     eff_comp = 1;
    T1 = 298;
    P1 = 101; %[kPa]
    R = 8.314; %[kJ/(kmol*K)]
    h1 = 0; %enthalpy of air = 0 at 25oC
    mdot = 6.0; %[kg/s]
    M_molar = 28.97; %[kg/kmol]

    %% Iteration
    error = [1000 1000 1000]; %Initialize error
    T2 = [298 0 1200]; %Initial guess range of T2
    tol = 0.2; %Tolerance of 0.2 oC

    i=1;
    while abs(error(2))>tol %|| (error(1)*error(3)<0)

        mid(i) = (T2(1)+T2(3))/2;
        T2(2) = mid(i);
        Tavg1 = (T1+T2)./2;

        %(ii-iv): Compute T2_calc1
        cp1 = 28.11 + Tavg1.*0.1967e-2 + Tavg1.^2.*0.4802e-5 + Tavg1.^3.*-1.966e-9;
        cv1 = cp1-R;
        k1 = cp1./cv1; %Array
        T2s_1 = T1*(P2/P1).^((k1-1)./k1);
        T2_calc1 = T1 + (T2s_1-T1)./eff_comp; %array

        %(v-vii): Compute T2_calc2, Recomputation using T2_calc1 for Tavg
        Tavg2 = (T1+T2_calc1)./2; %Array

        %(ii-iv): Compute T2_calc1
        cp2 = 28.11 + Tavg2.*0.1967e-2 + Tavg2.^2.*0.4802e-5 + Tavg2.^3.*-1.966e-9;
        cv2 = cp2-R;
        k2 = cp2./cv2; %Array
        T2s_2 = T1*(P2/P1).^((k2-1)./k2);
        T2_calc2 = T1 + (T2s_2-T1)./eff_comp; 

        %Error for Higher and Lower Bounds
        error = T2_calc2-T2_calc1; %array

        if error(1)*error(2)<0
            T2(3) = T2(2);
        else
            T2(1) = T2(2);
        end

        %For plotting/debugging
        T2_result1(i) = T2_calc1(2);
        T2_result2(i) = T2_calc2(2);
        errorPlot(i) = error(2);
        cp_result = cp2(2);
        i=i+1;
    end
    
    T2_out = T2_result2(end);
    h2_out = cp_result*(T2_out-T1)+h1; %[kJ/kmol]
    W_req = mdot*cp_result/M_molar*(T2_out-T1); %[kW]
    
    % Plots
%     plot(T2_result1); hold on; plot(T2_result2);
%     xlabel("Iterations i");
%     ylabel("Temperature [K]");
    
 end