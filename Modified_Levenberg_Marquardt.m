%This function performs the Modified Levenberg-Marquardt method of solving the
%inverse kinematics problem
%function Modified_Levenberg_Marquardt
clc;

disp('Modified Levenberg-Marquardt Method selected');
syms q1 q2 q3 q4 q5

%select the degree of freedom(DOF)
disp(['Select Manipulators DOF' newline '2. 2-DOF' newline '3. 3-DOF' newline '4. 4-DOF' newline '5. 5-DOF']);
dof = str2double(input('DOF: ','s'));

%Receive inputs
goal_point = input('Goal point (xf) = ');
xf = transpose(goal_point);
xi = 0.01; %ξ
error_margin = input('Error margin (ε) = ');
a1 = input('link length 1 = ');
a2 = input('link length 2 = ');

iteration_counter = 0;
%definitions(iteration_counter);
iterations = zeros(1,iteration_counter);
distance = zeros(1,iteration_counter);
x_position = zeros(1,iteration_counter);
y_position = zeros(1,iteration_counter);
kin = zeros(2,iteration_counter);
alpha_rad = zeros(1,iteration_counter);
alpha_grad = zeros(1,iteration_counter);
alpha_deg = zeros(1,iteration_counter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The section below is a switch case function to select which method to
%execute
switch dof
    case 2
        init_config = input('Initial Config. (q0) = ');
        q0 = transpose((init_config)*(pi/180));
        config = zeros(length(q0),iteration_counter);
        kinematics = [a1*cos(q1) + a2*cos(q1+q2); a1*sin(q1) + a2*sin(q1+q2)];
          
        Ja = jacobian(kinematics,[q1;q2]);  Jt = transpose(Ja);
        %M_matrix = Ja*Jt;
        
        %calculate qi
        q = q0;
        tic
        for i = 1:inf
            kinematics_ = subs(kinematics, [q1 q2], q(:,i).');
            Ja_ = subs(Ja, [q1 q2], q(:,i).');
            Jt_ = subs(Jt, [q1 q2], q(:,i).');
            I_matrix = eye(length(diag(inv(Ja_*Jt_))));
                        
            q(:,i+1) = q(:,i) + xi*Jt_*inv(diag(Ja_*Jt_) + I_matrix)*(xf - kinematics_);
            toc

            config(:,i) = q(:,i)*(180/pi);
            kin(:,i) = eval(kinematics_);
            
            iterations(i) = i;
            x_position(i) = eval(kinematics_(1));%saves value of x in kinematics_
            y_position(i) = eval(kinematics_(2));
            iteration_counter = iteration_counter + 1;
            
            dist = norm(xf - kinematics_);

            print_2DOF(i, dist, kinematics_, config)
            
            if dist < error_margin
                break ;
            end           
        end
        x_ = goal_point(:,1);
        y_ = goal_point(:,2);
        
        pose1 = goal_point(:,1) - x_position(1,1);
        pose2 = y_position(1,1) - y_position(1,:);
        pose3 = x_position(1,1) - x_position(1,:);
        pose4 = goal_point(:,2) - y_position(1,1);
        
        distance = abs((pose1*pose2)-(pose3*pose4))/(sqrt((pose1)^2)+(pose4)^2);
        average_distance = sum(distance)/iteration_counter;
        fprintf('Average Distance = %f\n',average_distance);

        w1 = zeros(2,iteration_counter);
        w2 = zeros(2,iteration_counter);
        w1w2 = zeros(1,iteration_counter);

        for n = 1:length(kin)-1
            w1(:,n) = kin(:,n+1) - kin(:,n);
            w2(:,n) = xf - kin(:,n);
            w1w2(:,:) = dot(w1,w2);
            w1_n = norm(w1);
            w2_n = norm(w2);
            alpha_deg = (acos((w1w2)/dot(w1_n,w2_n)))/(pi/180);
            alpha_rad = acos((w1w2)/dot(w1_n,w2_n));
            alpha_grad = alpha_rad*(pi/200);
        end
        calculate_print_angle_ave(iteration_counter,alpha_rad,alpha_deg,alpha_grad);
        
        figure (1)
        plot(iterations,x_position, 'b.-');
        xlim([0 (iteration_counter)]);
        ylim([(x_ -0.5) inf]);
        hold on
        plot(iteration_counter,x_,'Marker','*', 'Markersize',8)
        grid on
        xlabel('Number of Iterations')
        ylabel('x^{ith} Position')
        title('Modified Levenberg-Marquardt 2-DOF')
        plots(iterations,x_position,y_position,distance,alpha_deg,alpha_rad,alpha_grad,x_,y_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        a3 = input('link length 3 = ');
        init_config = input('Initial Config. (q0) = ');
        q0 = transpose((init_config)*(pi/180));
        config = zeros(length(q0),iteration_counter);
        kinematics = [a1*cos(q1) + a2*cos(q1+q2) + a3*cos(q1+q2+q3); a1*sin(q1) + a2*sin(q1+q2) + a3*sin(q1+q2+q3)];
                
        Ja = jacobian(kinematics,[q1;q2;q3]);  Jt = transpose(Ja);
        %M_matrix = Ja*Jt; 
        
        q = q0;
        tic
        for i = 1:inf
            kinematics_ = subs(kinematics, [q1 q2 q3], q(:,i).');
            Ja_ = subs(Ja, [q1 q2 q3], q(:,i).');
            Jt_ = subs(Jt, [q1 q2 q3], q(:,i).');
            I_matrix = eye(length(diag(inv(Ja_*Jt_))));

            q(:,i+1) = q(:,i) + xi*Jt_*inv(diag(Ja_*Jt_) + I_matrix)*(xf - kinematics_);
            toc

            config(:,i) = q(:,i)*(180/pi);
            kin(:,i) = eval(kinematics_);

            iterations(i) = i;
            x_position(i) = eval(kinematics_(1));%saves value of x in kinematics_
            y_position(i) = eval(kinematics_(2));
            iteration_counter = iteration_counter + 1;

            dist = norm(xf - kinematics_);

            print_3DOF(i, dist, kinematics_, config)
            
            if dist < error_margin
                break ;
            end            
        end
        x_ = goal_point(:,1);
        y_ = goal_point(:,2);
        
        pose1 = goal_point(:,1) - x_position(1,1);
        pose2 = y_position(1,1) - y_position(1,:);
        pose3 = x_position(1,1) - x_position(1,:);
        pose4 = goal_point(:,2) - y_position(1,1);
        
        distance = abs((pose1*pose2)-(pose3*pose4))/(sqrt((pose1)^2)+(pose4)^2);
        average_distance = sum(distance)/iteration_counter;
        fprintf('Average Distance = %f\n',average_distance);

        w1 = zeros(2,iteration_counter);
        w2 = zeros(2,iteration_counter);
        w1w2 = zeros(1,iteration_counter);

        for n = 1:length(kin)-1
            w1(:,n) = kin(:,n+1) - kin(:,n);
            w2(:,n) = xf - kin(:,n);
            w1w2(:,:) = dot(w1,w2);
            w1_n = norm(w1);
            w2_n = norm(w2);
            alpha_deg = (acos((w1w2)/dot(w1_n,w2_n)))/(pi/180);
            alpha_rad = acos((w1w2)/dot(w1_n,w2_n));
            alpha_grad = alpha_rad*(pi/200);
        end
        calculate_print_angle_ave(iteration_counter,alpha_rad,alpha_deg,alpha_grad);

        figure (1)
        plot(iterations,x_position, 'b.-');
        xlim([0 (iteration_counter)]);
        ylim([(x_ -0.5) inf]);
        hold on
        plot(iteration_counter,x_,'Marker','*', 'Markersize',8)
        grid on
        xlabel('Number of Iterations')
        ylabel('x^{ith} Position')
        title('Modified Levenberg-Marquardt 3-DOF')
        plots(iterations,x_position,y_position,distance,alpha_deg,alpha_rad,alpha_grad,x_,y_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        a3 = input('link length 3 = ');
        a4 = input('link length 4 = ');
        init_config = input('Initial Config. (q0) = ');
        q0 = transpose((init_config)*(pi/180));
        config = zeros(length(q0),iteration_counter);
        kinematics = [a1*cos(q1) + a2*cos(q1+q2) + a3*cos(q1+q2+q3) + a4*cos(q1+q2+q3+q4); ...;
                      a1*sin(q1) + a2*sin(q1+q2) + a3*sin(q1+q2+q3) + a4*sin(q1+q2+q3+q4)];
                
        Ja = jacobian(kinematics,[q1;q2;q3;q4]);  Jt = transpose(Ja);
        %M_matrix = Ja*Jt; 

        q = q0;
        tic
        for i = 1:inf
            kinematics_ = subs(kinematics, [q1 q2 q3 q4], q(:,i).');
            Ja_ = subs(Ja, [q1 q2 q3 q4], q(:,i).');
            Jt_ = subs(Jt, [q1 q2 q3 q4], q(:,i).');
            I_matrix = eye(length(diag(inv(Ja_*Jt_))));

            q(:,i+1) = q(:,i) + xi*Jt_*inv(diag(Ja_*Jt_) + I_matrix)*(xf - kinematics_);
            toc

            config(:,i) = q(:,i)*(180/pi);
            kin(:,i) = eval(kinematics_);

            iterations(i) = i;
            x_position(i) = eval(kinematics_(1));%saves value of x in kinematics_
            y_position(i) = eval(kinematics_(2));
            iteration_counter = iteration_counter + 1;

            dist = norm(xf - kinematics_);

            print_4DOF(i, dist, kinematics_, config)
            
            if dist < error_margin
                break ;
            end            
        end
        x_ = goal_point(:,1);
        y_ = goal_point(:,2);
        
        pose1 = goal_point(:,1) - x_position(1,1);
        pose2 = y_position(1,1) - y_position(1,:);
        pose3 = x_position(1,1) - x_position(1,:);
        pose4 = goal_point(:,2) - y_position(1,1);
        
        distance = abs((pose1*pose2)-(pose3*pose4))/(sqrt((pose1)^2)+(pose4)^2);
        average_distance = sum(distance)/iteration_counter;
        fprintf('Average Distance = %f\n',average_distance);

        w1 = zeros(2,iteration_counter);
        w2 = zeros(2,iteration_counter);
        w1w2 = zeros(1,iteration_counter);

        for n = 1:length(kin)-1
            w1(:,n) = kin(:,n+1) - kin(:,n);
            w2(:,n) = xf - kin(:,n);
            w1w2(:,:) = dot(w1,w2);
            w1_n = norm(w1);
            w2_n = norm(w2);
            alpha_deg = (acos((w1w2)/dot(w1_n,w2_n)))/(pi/180);
            alpha_rad = acos((w1w2)/dot(w1_n,w2_n));
            alpha_grad = alpha_rad*(pi/200);
        end
        calculate_print_angle_ave(iteration_counter,alpha_rad,alpha_deg,alpha_grad);
        
        figure (1)
        plot(iterations,x_position, 'b.-');
        xlim([0 (iteration_counter)]);
        ylim([(x_ -0.5) inf]);
        hold on
        plot(iteration_counter,x_,'Marker','*', 'Markersize',8)
        grid on
        xlabel('Number of Iterations')
        ylabel('x^{ith} Position')
        title('Modified Levenberg-Marquardt 4-DOF')
        plots(iterations,x_position,y_position,distance,alpha_deg,alpha_rad,alpha_grad,x_,y_)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        a3 = input('link length 3 = ');
        a4 = input('link length 4 = ');
        a5 = input('link length 5 = ');
        init_config = input('Initial Config. (q0) = ');
        q0 = transpose((init_config)*(pi/180));
        config = zeros(length(q0),iteration_counter);
        kinematics = [a1*cos(q1) + a2*cos(q1+q2) + a3*cos(q1+q2+q3) + a4*cos(q1+q2+q3+q4) + a5*cos(q1+q2+q3+q4+q5); ...;
                      a1*sin(q1) + a2*sin(q1+q2) + a3*sin(q1+q2+q3) + a4*sin(q1+q2+q3+q4) + a5*sin(q1+q2+q3+q4+q5)];
                
        Ja = jacobian(kinematics,[q1;q2;q3;q4;q5]);  Jt = transpose(Ja);
        %M_matrix = Ja*Jt; 
               
        q = q0;
        tic
        for i = 1:inf
            kinematics_ = subs(kinematics, [q1 q2 q3 q4 q5], q(:,i).');
            Ja_ = subs(Ja, [q1 q2 q3 q4 q5], q(:,i).');
            Jt_ = subs(Jt, [q1 q2 q3 q4 q5], q(:,i).');
            I_matrix = eye(length(diag(inv(Ja_*Jt_))));

            q(:,i+1) = q(:,i) + xi*Jt_*inv(diag(Ja_*Jt_) + I_matrix)*(xf - kinematics_);
            toc

            config(:,i) = q(:,i)*(180/pi);
            kin(:,i) = eval(kinematics_);

            iterations(i) = i;
            x_position(i) = eval(kinematics_(1));%saves value of x in kinematics_
            y_position(i) = eval(kinematics_(2));
            iteration_counter = iteration_counter + 1;

            dist = norm(xf - kinematics_);

            print_5DOF(i, dist, kinematics_, config)
            
            if dist < error_margin
                break ;
            end          
        end
        x_ = goal_point(:,1);
        y_ = goal_point(:,2);
        
        pose1 = goal_point(:,1) - x_position(1,1);
        pose2 = y_position(1,1) - y_position(1,:);
        pose3 = x_position(1,1) - x_position(1,:);
        pose4 = goal_point(:,2) - y_position(1,1);
        
        distance = abs((pose1*pose2)-(pose3*pose4))/(sqrt((pose1)^2)+(pose4)^2);
        average_distance = sum(distance)/iteration_counter;
        fprintf('Average Distance = %f\n',average_distance);

        w1 = zeros(2,iteration_counter);
        w2 = zeros(2,iteration_counter);
        w1w2 = zeros(1,iteration_counter);

        for n = 1:length(kin)-1
            w1(:,n) = kin(:,n+1) - kin(:,n);
            w2(:,n) = xf - kin(:,n);
            w1w2(:,:) = dot(w1,w2);
            w1_n = norm(w1);
            w2_n = norm(w2);
            alpha_deg = (acos((w1w2)/dot(w1_n,w2_n)))/(pi/180);
            alpha_rad = acos((w1w2)/dot(w1_n,w2_n));
            alpha_grad = alpha_rad*(pi/200);
        end
        calculate_print_angle_ave(iteration_counter,alpha_rad,alpha_deg,alpha_grad);
        
        figure (1)
        plot(iterations,x_position, 'b.-');
        xlim([0 (iteration_counter)]);
        ylim([(x_ -0.5) inf]);
        hold on
        plot(iteration_counter,x_,'Marker','*', 'Markersize',8)
        grid on
        xlabel('Number of Iterations')
        ylabel('x^{ith} Position')
        title('Modified Levenberg-Marquardt 5-DOF')
        plots(iterations,x_position,y_position,distance,alpha_deg,alpha_rad,alpha_grad,x_,y_)

    otherwise
        disp('Incorrect selection. Try again')
        Jacobian_Transpose_Method
end
