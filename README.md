# A_Comparison_of_Selected_Method_of_Inverse_Kinematics_for_Robot_Manipulators
This contains my Msc Thesis in Automation and Robotics. The project implemets three methods of solving the inverse Kinematics problem: Jacobian Transpose, Jacobian Pseudo-Inverse and the Modified Levenberg-Marquardt Methods.

The **run.m** file is executed first, then the method is selected.  Note that any value outside the displayed values are ignored and it prompts a re-entry of another value. After selecting the method, the DOF (Degree of Freedom is selected). In this case, we are limited to 2, 3, 4 and 5 DOF. To increase the DOF, the kinematics equation in the code should be modified appropriately, and other dependent parameters tweaked too. Other parameters as goal point, error margin, configuration, manipulator arm length are also inputed before the algorithm starts to run.

The Jacobian pseudo-inverse is the best and gives the most optimal solution while the Jacobian transpose is the fastest among the three methods.

To further understand how it works and how the algorithm was developed, go through the PDF document attached.
