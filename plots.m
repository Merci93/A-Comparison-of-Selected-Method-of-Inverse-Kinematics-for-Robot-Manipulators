function [] = plots(iterations,x_position,y_position,distance,alpha_deg,alpha_rad,alpha_grad,x_,y_,point2,point1)

        figure (2)
        plot(x_position,y_position, 'b.-');
        title('End effector position');
        xlabel('X-location');
        ylabel('Y-location');
        hold on
        plot(x_,y_,'Marker','*', 'Markersize',8)
        %hold on
        %plot(point2,point1, 'r.-')
        grid on

        figure (3)
        plot(iterations,distance, 'b.-')
        xlabel("Number of iterations")
        ylabel("Distance between point and straight line")
        title("Distance Plot")
        grid on

        figure(4)
        plot(iterations,alpha_grad, 'r.-');
        xlabel("Number of iterations");
        ylabel("Angle Quality Measure (α) in gradians (grad)");
        title("Angle Quality Measure Plot in Gradians");
        grid on

%         figure(5)
%         plot(iterations,alpha_deg, 'r.-');
%         xlabel("Number of iterations");
%         ylabel("Angle Quality Measure (α) in degrees (^o)");
%         title("Plot of Angle Quality Measure in Degrees");
%         grid on
% 
%         figure(6)
%         plot(iterations,alpha_rad, 'r.-');
%         xlabel("Number of iterations");
%         ylabel("Angle Quality Measure (α) in radians (rad)");
%         title("Plot of Angle Quality Measure in Radians");
%         grid on

end