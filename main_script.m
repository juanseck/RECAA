%  IEEE Access, DOI: http://
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj
% The initial parameters that you need are:
%__________________________________________
% SmartCells_no = number of smart-cells
% Neighbors_no = number of neighbors
% Max_iteration = maximum number of iterations
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% dim = number of variables to be tunned
% fobj = @YourCostFunction
% flag = print the optimization process every flag steps, print nothing if flag=0 
%
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers
%
% To run RECAA: [min_value,position_vector,convergence_curve]=RECAA(SmartCells_no,Neighbors_no,Max_iteration,lb,ub,dim,fobj,flag)
%______________________________________________________________________________________________


clear all
clc

Function_name='F3';     %Name of the test function that can be from F1 to F50  

SmartCells_no=12;       % Number of smart-cells
Neighbors_no=6;         % Number of neighbors for each smart-cell
Max_iteration=500;      % Maximum numbef of iterations
Solution_no=SmartCells_no*Neighbors_no;

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=benchmark_functions(Function_name);
%dim=500;

% Execute optimization algorithm
% Algoritmo MmCAA
[min_value,position_vector,convergence_curve]=RECAA(SmartCells_no,Neighbors_no,Max_iteration,lb,ub,dim,fobj,100);

figure(1)
if sum(convergence_curve<0) >0
    plot(1:Max_iteration,convergence_curve,'-s','Color','r','LineWidth',1.5,'MarkerSize',10,'MarkerIndices',1:50:Max_iteration)
else
    semilogy(1:Max_iteration,convergence_curve,'-s','Color','r','LineWidth',1.5,'MarkerSize',10,'MarkerIndices',1:50:Max_iteration)
end
legend('RECAA')
title(Function_name,'Fontsize',14)
xlabel('Iterations','Fontsize',13);
ylabel('Best fitness','Fontsize',13);
axis tight

display(['The best optimal value of the objective funciton found by RECAA is : ', num2str(min_value)]);

