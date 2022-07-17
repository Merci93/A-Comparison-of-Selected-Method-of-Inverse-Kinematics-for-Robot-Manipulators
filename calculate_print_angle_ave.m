function [] = calculate_print_angle_ave(iteration_counter,alpha_rad,alpha_deg,alpha_grad)

radians = sum(alpha_rad)/iteration_counter;
degrees = sum(alpha_deg)/iteration_counter;
gradians = sum(alpha_grad)/iteration_counter;

fprintf("Average angle in radians: %f\n", radians);
fprintf("Average angle in degrees: %f\n", degrees);
fprintf("Average angle in gradians: %f\n", gradians);

end
