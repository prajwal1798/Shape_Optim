function [Cl, Cd] = computeLiftDrag(uwallData, lwallData)
    % FARFIELD CONSTANTS
    alpha_deg = 0;                    % AOA IN DEGREES
    rho_inf = 1;                      % FREESTREAM RHO
    u_inf = 1;                        % FREESTREAM U
    q_inf = 0.5 * rho_inf * u_inf^2;  % DYN. PRESSURE
    S_ref = 1;                        % REF. AREA

    % FREESTREAM VECTORS
    n_inf = [-sind(alpha_deg), cosd(alpha_deg)];
    t_inf = [cosd(alpha_deg), sind(alpha_deg)];

    lwallData = flipud(lwallData); 
    fullData = [uwallData; lwallData];

    % EXTRACT DATA
    X = fullData(:, 1);     % X COORD
    Y = fullData(:, 2);     % X COORD
    p = fullData(:, 3);     % Pressure    
    tau1 = fullData(:, 4);  % TAU1
    tau2 = fullData(:, 5);  % TAU2
    n1 = fullData(:, 6);    % N1
    n2 = fullData(:, 7);    % N2 

    % INITIALIZE INTEGRALS
    integral_lift = 0;
    integral_drag = 0;

    % CALCULATE EDGES
    num_edges = length(p) - 1;
    for i = 1:num_edges
        
        x1 = X(i); y1 = Y(i);
        x2 = X(i + 1); y2 = Y(i + 1);
        p1 = p(i); p2 = p(i + 1);
        tau1_1 = tau1(i); tau1_2 = tau1(i + 1);
        tau2_1 = tau2(i); tau2_2 = tau2(i + 1);
        n1_1 = n1(i); n1_2 = n1(i + 1);
        n2_1 = n2(i); n2_2 = n2(i + 1);

        % COMPUTE JACOBIAN
        J = sqrt(0.25 * ((x2 - x1)^2 + (y2 - y1)^2));

        % STRESS TENSOR COMPONENTS
        sigma1_w1 = -p1 * n1_1 + tau1_1;
        sigma1_w2 = -p2 * n1_2 + tau1_2;
        sigma2_w1 = -p1 * n2_1 + tau2_1;
        sigma2_w2 = -p2 * n2_2 + tau2_2;

        % AVERAGING OVER EDGE
        sigma1_w = (sigma1_w1 + sigma1_w2) / 2;
        sigma2_w = (sigma2_w1 + sigma2_w2) / 2;

        % NUM. INTEGRAL (1 POINT GAUSS QUAD.) 
        w = 2;  % WEIGHT
        integral_lift = integral_lift + ((sigma1_w * n_inf(1) + sigma2_w * n_inf(2)) * w * J);
        integral_drag = integral_drag + ((sigma1_w * t_inf(1) + sigma2_w * t_inf(2)) * w * J);
    end
    
    % CD, CL
    Cl = -integral_lift / (q_inf * S_ref);
    Cd = -integral_drag / (q_inf * S_ref);
end
