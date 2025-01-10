function flippedControlPoints = flipControlPoints(controlPoints)
    % Reshape controlPoints to ensure it's a column vector of length 22
    controlPoints = reshape(controlPoints, 22, 1);

    % Preallocate flippedControlPoints as a column vector of 22 elements
    flippedControlPoints = zeros(22, 1);

    % Apply flipping for the upper surface (first 12 control points)
    flippedControlPoints(1) = controlPoints(11);  % Row 11 to 1
    flippedControlPoints(2) = controlPoints(12);  % Row 12 to 2
    flippedControlPoints(3) = controlPoints(9);   % Row 9 to 3
    flippedControlPoints(4) = controlPoints(10);  % Row 10 to 4
    flippedControlPoints(5) = controlPoints(7);   % Row 7 to 5
    flippedControlPoints(6) = controlPoints(8);   % Row 8 to 6
    flippedControlPoints(7) = controlPoints(5);   % Row 5 to 7
    flippedControlPoints(8) = controlPoints(6);   % Row 6 to 8
    flippedControlPoints(9) = controlPoints(3);   % Row 3 to 9
    flippedControlPoints(10) = controlPoints(4);  % Row 4 to 10
    flippedControlPoints(11) = controlPoints(1);  % Row 1 to 11
    flippedControlPoints(12) = controlPoints(2);  % Row 2 to 12

    % Apply flipping for the lower surface (last 10 control points)
    flippedControlPoints(13) = controlPoints(21); % Row 21 to 13
    flippedControlPoints(14) = controlPoints(22); % Row 22 to 14
    flippedControlPoints(15) = controlPoints(19); % Row 19 to 15
    flippedControlPoints(16) = controlPoints(20); % Row 20 to 16
    flippedControlPoints(17) = controlPoints(17); % Row 17 remains 17
    flippedControlPoints(18) = controlPoints(18); % Row 18 remains 18
    flippedControlPoints(19) = controlPoints(15); % Row 15 to 19
    flippedControlPoints(20) = controlPoints(16); % Row 16 to 20
    flippedControlPoints(21) = controlPoints(13); % Row 13 to 21
    flippedControlPoints(22) = controlPoints(14); % Row 14 to 22

    % Reshape to (1, 22) for NN compatibility
    flippedControlPoints = reshape(flippedControlPoints, 1, 22);
end
