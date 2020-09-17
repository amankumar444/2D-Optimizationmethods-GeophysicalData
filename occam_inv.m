clear all;
clc;

disp('Determine number of layers:  from 6 - 30:  ');
n = input('Enter total number of layers:: ');
% To determine if we want to see the plot or not
disp('Want Plot of every Model:: 1 - Yes or 2 - No:  ');
fig_input = input('1::Yes, 2::No: ');
% To determine if we want to save the results of each layer or not
disp('Want .dat file of every Model:: 1 - Yes or 2 - No');
ok = input('1::Yes, 2::No: ');

n_d = input('Enter total number of soundings to check:: ');

str_var = input('Character of datafile before numeric:: ','s');
%**************************************************************************
rho=100*ones(1,n);%30 layer initial model
Geom=(1)*[1:length(rho)-1];
%**************************************************************************
var = 4;
d=var*[Geom];
m_finforplot = [];
Resistivity = rho;
i = sqrt(-1);
RMS = [];
Td_res = [];
Td_total = [];
for y = 1:1:n_d
    m=[rho]';%occam
    d0 = [];
    y
    str2fu = {'DAta'};
    [path] = path_check(pwd,str2fu);
    [f,re,Im,r,h] = read_data(path,strcat(str_var,num2str(y),'_marq.dat'));

    d0 = [d0                                               % initialization of d0
         re
         Im];

    dp=forward_HEM(rho,d,h,r,f);
    W=zeros(length(d0));
    for i=1:length(d0)
    W(i,i)=50;                                              % 1/(0.02*dobs(i));
    end
    u = 0;                                                  % mew value in the equation
    P = delta(length(m));                                   % first order roughness
    n = 1;
    while(n<6)
        n = n+1;
        [Resistivity,J,dp]=Jack(m,d,h,r,f,m);               % jacobian
        d2 = d0-dp+J*m;
        temp1 = u*P'*P;
        temp2 = (W*J)'*W*J;
        mnew = pinv(temp1+temp2)*(W*J)'*W*d2;
        rms = norm(d0-dp)/norm(d0);
        m = mnew;       
    end

    RMS(y) = norm(dp-d0)/norm(d0);
    Robs=d0(1:length(d0)/2);
    Rp=dp(1:length(d0)/2);
    Imobs=d0(length(d0)/2+1:length(d0));
    Imp=dp(length(dp)/2+1:length(dp));
    
    m_finforplot = [m_finforplot
        m];
    
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
        str2fu = {'Res','datafile_res','occam'};
        [path] = path_check(pwd,str2fu);
        name = strcat('Occam_',num2str(y),'.dat');
        fileName = fopen(fullfile(path, name),'w');
        % Writing data to file
        fprintf(fileName, '%f %f %f %f %f\n',fin_res);
        % Closing
        fclose(fileName); 
    end
    %**********************************************************************
    %conversion to apparent resistivity
    
    
    % Plots  of Predicted real and Imaginary plot
        
        F_depth = [];
        F_Resistivity = [];
        
        F_depth = [0];
        F_Resistivity = [Resistivity(1)];
        Td_res = [Td_res
            F_Resistivity];
        
        
        t_change = 1;
        count = 1;
        for i = 1:1:var*length(Resistivity)
           F_depth = [F_depth
               -1.*i];
           if t_change ~= var
               F_Resistivity = [F_Resistivity
                   Resistivity(count)];
               t_change = t_change + 1;
           else
               t_change = 1;
               F_Resistivity = [F_Resistivity
                   Resistivity(count)];
               count = count+1;
           end        
        end
        Td_total = [Td_total
        F_Resistivity(2:end)'];
    
        figure('Visible', 'off');
        
        if(fig_input == 1)
        
        ax1 = subplot(15,11,[23,34,45,56,2,3,4,5,24,25,26,27,35,36,37,38,46,47,48,49,57,58,59,60]);
        loglog(f,Robs,f,Rp,'-ro');
        xlabel('Frequency');
        ylabel('Inphase');
        title('Occam Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax2 = subplot(15,11,[100,111,122,133,144,155,101,102,103,104,112,113,114,115,123,124,125,126,134,135,136,137,145,146,147,148,156,157,158,159]);
        loglog(f,Imobs,f,Imp,'-mo');
        xlabel('Frequency');
        ylabel('quadrature');
        title('Occam Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax3 = subplot(15,11,[18,19,20,21,22,29,30,31,32,40,41,42,43,51,52,53,54,62,63,64,65,73,74,75,76,84,85,86,87,95,96,97,98,106,107,108,109,117,118,119,120,128,129,130,131,139,140,141,142,150,151,152,153,161,162,163,164]);
        plot(F_Resistivity,F_depth,'r');
        LineWidth = 1.4;
        xlabel('Resistivity');
        ylabel('Depth');
        xlim([min(F_Resistivity)-10 max(F_Resistivity)+10]);
        title({strcat('(RMSE =  ',num2str(RMS(y)),')'),'Resistivity vs Depth'});
        
        str2fu = {'Res','Plot','Occam'};
        [path] = path_check(pwd,str2fu);
        name = strcat('Occam_',num2str(y),'.jpg');
        saveas(gcf, fullfile(path, name), 'jpg');
        %saveas(gcf,);
        %%%code to save plots in a drive
        
    end
end
disp('Plots and .dat file for every model is saved in Res Folder.');
%**************************************************************************
%Final plot resistivity vs depth for determined number of soundings

figure('visible', 'on');
imagesc(Td_total');
colormap(parula);
colorbar;
str2fu = {'Res','Plot','Occam'};
[path] = path_check(pwd,str2fu);
name = strcat('Occam_',num2str(y+1),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');
%**************************************************************************
%RMSE plot

figure('visible', 'on');
plot(RMS,'-kx');
xlabel('Soundings');
ylabel('Root Mean Error');
title('RMSE Values Marq');

str2fu = {'Res','Plot','Occam'};
[path] = path_check(pwd,str2fu);
name = strcat('Occam_',num2str(y+2),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');