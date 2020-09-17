clear all 
clc


disp('Determine the cross rate:: in 0-1 range ');
crossrate = input('Enter crossrate btw 0-1:: ');

% To determine if we want to see the plot or not
disp('Want Plot of every Model:: 1 - Yes or 2 - No');
fig_input = input('1::Yes, 2::No: ');

% To determine if we want to save the results of each layer or not
disp('Want .dat file of every Model:: 1 - Yes or 2 - No');
ok = input('1::Yes, 2::No: ');
%***************************************************************************************
RMS = [];
depth_final = [];
Resistivity_final = [];
str_var = input('Enter character of datafile before numeric:: ','s');
n_d = input('Enter total number of soundings to check:: ');
m_finforplot = [];
q = 0;
Td_total = [];
for y = 1:n_d
    y
    d0 = [];
    popsize = 1000;                     %Population size is initialized
    generations = 10;                   %Total number of iterations
    str2fu = {'DAta'};                  %Selected path for accessing Data files
    [path] = path_check(pwd,str2fu);    % Merger of path of subdirectory path and overall path of main directory
    [f,re,Im,r,h] = read_data(path,strcat(str_var,num2str(y),'_marq.dat'));
    d0 = [d0
         re
         Im];
    %%crossrate=0.9 %Initial Population
    check = d0';
    pop = rand(popsize,5);
    pop(:,1)= (pop(:,1)*40 +30);
    pop(:,2)= (pop(:,2)*50 +100);
    pop(:,3)= (pop(:,3)*30 +30);
    pop(:,4)= (pop(:,4)*70+190);
    pop(:,5)= (pop(:,5)*70 +30);
    out = zeros(popsize,10);
    fit = zeros(popsize,1);
    
    %Genetic Algorithm optimization::
    
    fit=fitness(pop,popsize,check);
    [bestfit best]=max(fit);
    best=pop(best,:);
    
    % Model Application::
    for(k= 1:generations)
        s = select(fit);
        newpop=cross(pop,s,crossrate);
        fit=fitness(newpop,popsize,check);
        if (max(fit)>bestfit)
            [bestfit best]=max(fit);
             best=newpop(best,:);
        end
        pop=newpop;
        bestfit;
        best;
        
    end
    dp_min = [];
    depth_final_M = [];
    Resistivity_final_M = [];
    min_RMS=999999;
    for i = 1:popsize
        d = pop(i,1:2);
        rho= pop(i,3:5);
        dp = forward_HEM(rho,d,h,r,f);
        RMS_var = sum(((check'-dp)./(check')).^2);
        if(min_RMS>RMS_var)
            dp_min = dp;
            min_RMS = RMS_var;
            depth_final_M = d;
            Resistivity_final_M = rho;
        end
    end
    
    RMS(y) = min_RMS;
    depth_final = [depth_final
         depth_final_M];
    Resistivity_final = [Resistivity_final
        Resistivity_final_M];
    
    k = [Resistivity_final_M depth_final_M];

    if length(q == 0)
        m_finforplot = [m_finforplot
        k];
        q = length(k);
    else
        m_finforplot = [m_finforplot
        k];
    end
    
    
%****************************************************************************************
    Robs=d0(1:5);
    Rp=dp(1:5);
    Imobs=d0(6:10);
    Imp=dp(6:10);
%****************************************************************************************
    %for saving .dat file
    if ok == 1
    % Creating matrix to save results in .dat extension::
        fin_res=[f' 
            Robs' 
            Rp' 
            Imobs' 
            Imp'];
        %%Exporting data to .dat
        % Creating file to be written
        str2fu = {'Res','datafile_res','GA'};
        [path] = path_check(pwd,str2fu);
        name = strcat('GA_',num2str(crossrate),'_',num2str(y),'.dat');
        fileName = fopen(fullfile(path, name),'w');
        % Writing data to file
        fprintf(fileName, '%f %f %f %f %f\n',fin_res);
        % Closing
        fclose(fileName); 
    end

    %%% Resistivity versus Depth Plot of the BEST RMSE from 1000 iteration for each Sounding
    
        f_depth = depth_final_M;
        f_rho = Resistivity_final_M;
        F_depth = [0];
        F_Resistivity = [f_rho(1)];
        
%********************************************************************************************
        t_change = 1;
        k = f_depth(t_change);
        for i = 1:1:300
            if k > abs(F_depth(i))
                F_depth = [F_depth
                    -1.*i];
                F_Resistivity = [F_Resistivity
                    f_rho(t_change)];
            else
                t_change = t_change+1;
                if t_change > length(f_depth)
                    k = k*1000;
                else
                    k = f_depth(t_change);
                end
                F_depth = [F_depth
                    -1.*i];
                F_Resistivity = [F_Resistivity
                    f_rho(t_change)];
            end
        end
        
        Td_total = [Td_total
            F_Resistivity(2:end)'];
        
        % Graph Plot::
        if fig_input == 1
        figure('visible', 'off');
        ax1 = subplot(15,11,[23,34,45,56,2,3,4,5,24,25,26,27,35,36,37,38,46,47,48,49,57,58,59,60]);
        loglog(f,Robs,f,Rp,'-ro');
        xlabel('Frequency');
        ylabel('Inphase');
        title('GA Inversion');
        legend('Observed data','Predicted data','FontSize',2,'Location','south');
        ax2 = subplot(15,11,[100,111,122,133,144,155,101,102,103,104,112,113,114,115,123,124,125,126,134,135,136,137,145,146,147,148,156,157,158,159]);
        loglog(f,Imobs,f,Imp,'-mo');
        xlabel('Frequency');
        ylabel('quadrature');
        title('GA Inversion');
        legend('Observed data','Predicted data','FontSize',2,'Location','south');
        ax3 = subplot(15,11,[18,19,20,21,22,29,30,31,32,40,41,42,43,51,52,53,54,62,63,64,65,73,74,75,76,84,85,86,87,95,96,97,98,106,107,108,109,117,118,119,120,128,129,130,131,139,140,141,142,150,151,152,153,161,162,163,164]);
        plot(F_Resistivity,F_depth,'r');
        xlabel('Resistivity');
        ylabel('Depth');
        xlim([min(F_Resistivity)-10 max(F_Resistivity)+10]);
        LineWidth = 1.4;
        title({strcat('(RMSE =  ',num2str(min_RMS),')'),'Resistivity vs Depth'});
        
        str2fu = {'Res','Plot','GA'};
        [path] = path_check(pwd,str2fu);
        name = strcat('GA_',num2str(crossrate),'_',num2str(y),'.jpg');
        saveas(gcf, fullfile(path, name), 'jpg');
    end
    
end
%**************************************************************************
[M,I] = min(RMS);
disp('Best Root Mean Error');
M 
% disp('Best fit Depth');
depth_final(I,:)
disp('Best fit Resistivity');
Resistivity_final(I,:)

depth_final = depth_final(I,:);
Resistivity_final = Resistivity_final(I,:);

disp('Plots and .dat file for every model is saved in Res Folder.');
% %Final Resistivity plot corresponding to BEST RMSE score ::
% %************************************************************************
figure('visible', 'on');
imagesc(Td_total');
colormap(parula);
colorbar;
str2fu = {'Res','Plot','GA'};
[path] = path_check(pwd,str2fu);
name = strcat('GA_',num2str(y+3),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');

%**************************************************************************
%RMSE plot for every sounding

figure('visible','on');
plot(RMS,'-kx');
xlabel('Soundings');
ylabel('Root Mean Error');
title('RMSE Values GA');

str2fu = {'Res','Plot','GA'};
[path] = path_check(pwd,str2fu);
name = strcat('GA_',num2str(crossrate),'_',num2str(y+2),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');