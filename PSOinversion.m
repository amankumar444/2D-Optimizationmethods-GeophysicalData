clear all;
clc;

disp('determine number of layers:  from 3-5:  ');
n = input('Enter number of Layers: ');

% To determine if we want to see the plot or not
disp('Want Plot of every Model:: 1 - Yes or 2 - No');
fig_input = input('1::Yes, 2::No: ');

% To determine if we want to save the results of each layer or not
disp('Want .dat file of every Model:: 1 - Yes or 2 - No');
ok = input('1::Yes, 2::No: ');


n_d = input('Enter total number of soundings to check:: ');
str_var = input('Character of datafile before numeric:: ','s');
%********************************************************************
RMS = [];
depth_final=[];
Resistivity_final =[];
RMS_var = 0;
Min_RMS = 99999;
m_finforplot = [];
q=0;
Td_total = [];
for y = 1:n_d
    y
    str2fu = {'DAta'};
    [path] = path_check(pwd,str2fu);
    [f,re,Im,r,h] = read_data(path,strcat(str_var,num2str(y),'_marq.dat'));
    rho = 60.* ones(length(60),n);
    d = 50:50:(n-1)*50;
    %d=[1 thick];
    m=[rho d]';
    i=sqrt(-1);
    d0= [];
    d0 = [re
          Im];
       
    dp=forward_HEM(rho,d,h,r,f);
    m=PSOnew(m,d0,rho,h,r,f,d);
    dp=forward_HEM(m(1:length(rho)),m(length(rho)+1:length(m)),h,r,f);
    t = sqrt(norm((d0-dp)./d0)^2/length(d0));

    RMS_var = norm(dp-d0)/norm(d0);
    RMS(y) = RMS_var;
    Robs=d0(1:5);
    Rp=dp(1:5);
    Imobs=d0(6:10);
    Imp=dp(6:10);
    
    if length(q == 0)
            m_finforplot = [m_finforplot
                m];
            q = length(m);
    else
            m_finforplot = [m_finforplot
                m];
    end
    if(Min_RMS>RMS_var)
        Min_RMS = RMS_var;
        depth_final = m(length(rho)+1:length(m));
        Resistivity_final = m(1:length(rho));
    end  
    
    %for saving .dat file
    if ok == 1
    % Creating matrix to save results in .dat extension::
        fin_res=[f' 
            Robs' 
            Rp' 
            Imobs' 
            Imp'];
        %%Exporting data to .dat
        % Creating file to be written to
        str2fu = {'Res','datafile_res','PSO'};
        [path] = path_check(pwd,str2fu);
        name = strcat('PSO_',num2str(y),'.dat');
        fileName = fopen(fullfile(path, name),'w');
       
        % Writing data to file
        fprintf(fileName, '%f %f %f %f %f\n',fin_res);
        % Closing
        fclose(fileName); 
    end
    % Plots  of Predicted real and Imaginary plot
    
    f_depth = m(length(rho)+1:length(m));
    f_rho = m(1:length(rho));
        
    F_depth = [0];
    F_Resistivity = [f_rho(1)];
        
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
    
    if(fig_input == 1)
        figure('Visible', 'off');
        
        ax1 = subplot(15,11,[23,34,45,56,2,3,4,5,24,25,26,27,35,36,37,38,46,47,48,49,57,58,59,60]);
        loglog(f,Robs,f,Rp,'-ro');
        xlabel('Frequency');
        ylabel('Inphase');
        title('PSO Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax2 = subplot(15,11,[100,111,122,133,144,155,101,102,103,104,112,113,114,115,123,124,125,126,134,135,136,137,145,146,147,148,156,157,158,159]);
        loglog(f,Imobs,f,Imp,'-mo');
        xlabel('Frequency');
        ylabel('quadrature');
        title('PSO Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax3 = subplot(15,11,[18,19,20,21,22,29,30,31,32,40,41,42,43,51,52,53,54,62,63,64,65,73,74,75,76,84,85,86,87,95,96,97,98,106,107,108,109,117,118,119,120,128,129,130,131,139,140,141,142,150,151,152,153,161,162,163,164]);
        plot(F_Resistivity,F_depth,'r');
        LineWidth = 1.4;
        xlabel('Resistivity');
        ylabel('Depth');
        xlim([min(F_Resistivity)-10 max(F_Resistivity)+10]);
        title({strcat('(RMSE =  ',num2str(RMS(y)),')'),'Resistivity vs Depth'});
        
        str2fu = {'Res','Plot','PSo'};
        [path] = path_check(pwd,str2fu);
        name = strcat('PSO_',num2str(y),'.jpg');
        saveas(gcf, fullfile(path, name), 'jpg');

    end
end

%Final Resistivity plot::
%**************************************************************************
%%2D plot of all the soundings in a single heatplot::
figure('visible', 'on');
imagesc(Td_total');
colormap(parula);
colorbar;
str2fu = {'Res','Plot','PSo'};
[path] = path_check(pwd,str2fu);
name = strcat('PSO_',num2str(y+3),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');
%**************************************************************************


disp('Best Root Mean Error');
Min_RMS
% disp('Best fit Depth');
depth_final
disp('Best fit Resistivity');
Resistivity_final

disp('Plots and .dat file for every model is saved in Res Folder');
%**************************************************************************
%Code for plotting the results ::

figure('visible', 'on');
plot(RMS,'-kx');
xlabel('Soundings');
ylabel('Root Mean Error');
title('RMSE Values Marq');

str2fu = {'Res','Plot','PSo'};
[path] = path_check(pwd,str2fu);
name = strcat('PSO_',num2str(y+2),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');