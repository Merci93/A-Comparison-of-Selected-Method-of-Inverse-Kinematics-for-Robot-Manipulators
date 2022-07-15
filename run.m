% RUN: This script selects the program to be run and carries out the
% operation needed using the selected method.

clear global; %clear the global workspace
%clc; %clears the command window

addpath(pwd); %adds the matlab path
    
%This section requests for the function and system to be executed
disp(['Select method' newline '1. Jacobian Pseudo-Inverse Method' newline '2. Jacobian Transpose Method']);
disp(['3. Modified Levenberg-Marquardt' newline]);
selection = str2double(input('Select Method: ','s'));

%The section below is a switch case function to select which method to
%execute
switch selection
    case 1
        Jacobian_Pseudo_Inverse_Method
    case 2
        Jacobian_Transpose_Method
    case 3
        Modified_Levenberg_Marquardt
    otherwise
        fprintf('Incorrect selection. Try again\n\n')
        run
end
