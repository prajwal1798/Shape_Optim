function ub = getUpperBounds()
    % Upper bounds
    ub =  [0.4326; 0.0100; 0.1806; 0.0404; -0.0425; 0.0571; 
      -0.2484; 0.0627; -0.4397; 0.0474; -0.5000; 0.0175;
      0.4326; -0.0100; 0.1806; -0.0404; -0.0425; -0.0571;
      -0.2484; -0.0627; -0.4397; -0.0474] +  [
        5.3e-3; 2.7e-3; 1.5e-2; 7.4e-3; 2.5e-2; 1.2e-2;
        3.1e-2; 1.5e-2; 4.3e-2; 2.2e-2; 0; 2.2e-2;
        1.9e-2; 9.6e-3; 3.2e-2; 1.6e-2; 3.9e-5; 1.9e-5;
        4.5e-2; 2.3e-2; 4.3e-2; 2.2e-2
    ];
end
