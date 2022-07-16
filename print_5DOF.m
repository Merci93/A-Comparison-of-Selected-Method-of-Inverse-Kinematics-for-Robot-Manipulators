function [] = print_5DOF(i, dist, kinematics_, config)

fprintf('Number of iterations = %d\n',i);
fprintf('Distance between current position and goal = %f\n',dist)
fprintf('kinematics = [%f , %f]\n', kinematics_(1), kinematics_(2));
fprintf('Configurations (q) = [%f %f %f %f %f]\n\n',config(1,i),config(2,i),config(3,i),config(4,i),config(5,i));

end