%% Task 3 Problem Statement
% Compute cp_prod for product gas mixture as mole-fraction weighed average
% of the cp values. Assume stoichiometric conditions.

function cp_prod = Task3(a)
    %% Initializing Constants/Assumptions
    g = 0.25;
    Tavg = 600; %[K]

    % Test Case:
    % a_stoich = 4.76*(2+3*g);
    % a = a_stoich;

    mol_total = (3.76*a/4.76) + (a/4.76 - 3*g - 2) + (2 + 2*g) + (1 + 2*g);

    %% Molar Relations
    y_N2 = [3.76*a/4.76] / mol_total;
    y_O2 = [a/4.76 - 3*g - 2] / mol_total;
    y_H2O = [2 + 2*g] / mol_total;
    y_CO2 = [1 + 2*g] / mol_total;

    %   1) O2 
    %   2) H2O 
    %   3) CO2 
    %   4) N2

    a = [25.48 32.24 22.26 28.90];
    b = [1.520 0.1923 5.981 -0.1571].*10^-2;
    c = [-0.7155 1.055 -3.501 0.8081].*10^-5;
    d = [1.312 -3.595 7.469 -2.873].*10^-9;

    cp_O2 = a(1) + b(1)*Tavg + c(1)*Tavg^2 + d(1)*Tavg^3; %[kJ/(kmol*K)]
    cp_H2O = a(2) + b(2)*Tavg + c(2)*Tavg^2 + d(2)*Tavg^3;
    cp_CO2 = a(3) + b(3)*Tavg + c(3)*Tavg^2 + d(3)*Tavg^3;
    cp_N2 = a(4) + b(4)*Tavg + c(4)*Tavg^2 + d(4)*Tavg^3;

    cp_prod = y_H2O*cp_H2O + y_CO2*cp_CO2 + y_N2*cp_N2 + y_O2*cp_O2; %[kJ/(kmol*K)]

    %% Comparison to cp_air [Table pg. 7]
    cp_air = 28.11 + 0.1967e-2*Tavg + 0.4802e-5*Tavg^2 - 1.966e-9*Tavg^3;

    fprintf("At T = 600 K: \n   cp_prod = " + cp_prod + " [kJ/(kmol*K)] \n" + ...
        "   cp_air = " + cp_air + " [kJ/(kmol*K)]");
end
