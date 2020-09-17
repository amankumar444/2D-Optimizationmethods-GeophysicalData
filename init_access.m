%% NOTE: 
% 1. Please Go by the instructions provided in command Window once the code 
%    is running.
% 2. RMSE value has been calculated and plotted for each sounding.
% % 3. The Final Depth Versus Resistivity plot is calculated for the Least
%    RMSE score.
% 4. The .dat extension is saving Observed and Predicted Inphase and
%    Quadrature results for the corresponding Frequencies. 
% 5. Plot for each sounding gives Observed and Predicted results of  
%    Inphase and Quadrature values for differnt frequencies.
%% Marquadt Model:: 
% 1. Depth versus Resistivity plot for the sounding with Least RMSE score
%    is Plotted.
% 2. It is giving 3 layer model for the provided Dataset.
% 3. Files used:: "forward_HEM.m, marq_inv.m, Jack.m,RTE.m,read_data_m.m".
%% Occam Model:: 
% 1. The depth for each layer is predetermined as 2(two).
% 2. Files used::"delta.m,forward_HEM.m,Jack.m,RTE.m,occam_inv.m,read_data_m.m".
%% PSO::
% 1. There is problem in the code of PSO inversion model (it is giving
%    incorrect results of Depth).
% 2. Code is Running, but results are unsatisfactory.
% Files used:: "PSOinversion.m, PSOnew.m, read_data_m.m, RTE.m".
%% Genetic Algorithm::
% 1. It is also running for all 60 sounding Data and the resistivity vs depth graph is plotted for every sounding.
% 2. This Algorithm is iterating 1000 times in each sounding and the model with best RMSE score is selected from every 1000 iteration.
% 3. RMSE is Calculated for each iteration in each and every sounding
% 4. Best RMSE score for each sounding plotted.
% 5. Resistivity Versus Depth plot is done for sounding with least RMSE score 
%    outcome.
% 6. Model has considerably long running time
% 7. files used::
%          "cross.m,fitness.m,forward_HEM.m,GA_inversion.m,RTE.m,select.m, read_data_m.m "
%% Laterally Constrain::
% 1. 2 dimensional results has been calculated for this model.
% 2. The results for each sounding is saved.
% 3. files used::"forward_HEM.m,Jacob.m,Jacobian.m,read_data_m.m"
%%
%********************************************************************************************************
%% Final Code::
disp("Please select numbers between 1-5 for getting the result of different optimisation model");
disp("1:: Marquadt Inversion");
disp("2:: Occam Inversion");
disp("3:: PSO Inversion");
disp("4:: Genetic Algorithm Inversion");
disp("5:: Lateral constrain Inversion");
n = input("Enter the number between 1-5:: ");

switch n
    case 1
        marq_inv
    case 2
        occam_inv
    case 3
        PSOinversion
    case 4
        GA_inversion
    case 5
        norminvfinal
    otherwise
        disp('Error: X :');
end