disp('Determine number of layers:  from 3-5:  ');
l = input('Enter number of Layers: ');
disp('Determine initial density (rho):  ');
rho0 = input('Rho: ');
disp('Determine initial thickness(depth):  ');
d01 = input('depth in multiple of 10: ');

% To determine if we want to see the plot or not
disp('Want Plot of every Model:: 1 - Yes or 2 - No');
fig_input = input('1::Yes, 2::No: ');

% To determine if we want to save the results of each layer or not
disp('Want .dat file of every Model:: 1 - Yes or 2 - No');
ok = input('1::Yes, 2::No: ');
n = l;

n_d = input('Enter total number of soundings to check:: ');
disp('Enter data file initial name ');
str_var = input('Character of datafile before numeric:: ','s');
%********************************************************************
rho = rho0.* ones(length(60),l);
d=d01:d01:(n-1)*d01;
%********************************************************************
%Initialization of values:

m=[rho d]';
n_m = 2*l-1;
m_finforplot = [];
RMS = [];
Td_depth = [] ;
Td_Res = [];

for y = 1:1:n_d
    m=[rho d]';
    d0 = [];
    y
    str2fu = {'DAta'};
    [path] = path_check(pwd,str2fu);
    [f,re,Im,r,h] = read_data(path,strcat(str_var,num2str(y),'_marq.dat'));
    
    d0 = [d0
          re
          Im];
    dp=forward_HEM(rho,d,h,r,f);
    n = 0;
    u = 0.1;
    for i=1:length(d0)
        W(i,i)=1/d0(i);
    end
    m_pre = m;
    mini =m;
    tt =1;
    tt_pre = 2;
    tt_prep = 3;
    tt_min =4;
    %Optimization work::
    while(n<10&tt>0.065)
        
        n = n+1;
        jac = Jacobian(m(1:length(rho)),m(length(rho)+1:length(m)),h,r,f,m);
        q = inv(jac'*jac+u*eye(length(m)))*(jac')*(d0-dp);
        
        m = m + q;
        tt_prep = tt_pre;
        tt_pre = tt;
        tt = sqrt((norm((d0-dp)./d0)^2)/length(d0));
        if(tt<tt_min)
            tt_min = tt;
            m_minmisfit = m;
        end
        dp=forward_HEM(m(1:length(rho)),m(length(rho)+1:length(m)),h,r,f);
        rc = (tt-tt_pre)/(tt_prep-tt_pre);
        if(abs(rc)<0.01)
            u = 50;
        elseif(abs(rc)>0.25)
            u = 10;
        else
            u = 0.005;
        end
    end
    tt_min;
    m_finforplot = [m_finforplot
                    m_minmisfit];
                
    RMS(y) = norm(dp-d0)/norm(d0);
    dp = forward_HEM(m_minmisfit(1:length(rho)),m_minmisfit(length(rho)+1:length(m_minmisfit)),h,r,f);
    Robs=d0(1:5);
    Rp=dp(1:5);
    Imobs=d0(6:10);
    Imp=dp(6:10);
    
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
        str2fu = {'Res','datafile_res','laterally_constrain'};
        [path] = path_check(pwd,str2fu);
        name = strcat('lC_',num2str(y),'.dat');
        fileName = fopen(fullfile(path, name),'w');
        % Writing data to file
        fprintf(fileName,'%f %f %f %f %f\n',fin_res);
        % Closing
        fclose(fileName); 
    end
    %***************************************************************************
    % Conversion to apparent resistivity
    % rho_o=1./(2*pi*f*4*pi*10^-7).*((Robs+sqrt(-1)*Imobs).*(Robs-sqrt(-1)*Imobs));
    % rho_p=1./(2*pi*f*4*pi*10^-7).*((Rp+sqrt(-1)*Imp).*(Rp-sqrt(-1)*Imp));
    
    
        f_depth = m_minmisfit(length(rho)+1:length(m_minmisfit));
        f_rho = m_minmisfit(1:length(rho));
        
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
        
        Td_Res = [Td_Res
        F_Resistivity];
    
        if(fig_input == 1)
        figure('Visible', 'off');
        
        ax1 = subplot(15,11,[23,34,45,56,2,3,4,5,24,25,26,27,35,36,37,38,46,47,48,49,57,58,59,60]);
        loglog(f,Robs,f,Rp,'-ro');
        xlabel('Frequency');
        ylabel('Inphase');
        title('IC Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax2 = subplot(15,11,[100,111,122,133,144,155,101,102,103,104,112,113,114,115,123,124,125,126,134,135,136,137,145,146,147,148,156,157,158,159]);
        loglog(f,Imobs,f,Imp,'-mo');
        xlabel('Frequency');
        ylabel('quadrature');
        title('IC Inversion');
        legend('Observed data','Predicted data','FontSize',2);
        ax3 = subplot(15,11,[18,19,20,21,22,29,30,31,32,40,41,42,43,51,52,53,54,62,63,64,65,73,74,75,76,84,85,86,87,95,96,97,98,106,107,108,109,117,118,119,120,128,129,130,131,139,140,141,142,150,151,152,153,161,162,163,164]);
        plot(F_Resistivity,F_depth,'r');
        LineWidth = 1.4;
        xlabel('Resistivity');
        ylabel('Depth');
        xlim([min(F_Resistivity)-10 max(F_Resistivity)+10]);
        title({strcat('(RMSE =  ',num2str(RMS(y)),')'),'Resistivity vs Depth'});
        
        str2fu = {'Res','Plot','Laterally_constrain'};
        [path] = path_check(pwd,str2fu);
        name = strcat('lC_',num2str(y),'.jpg');
        saveas(gcf, fullfile(path, name), 'jpg');
        %saveas(gcf,);
        %%%code to save plots in a drive
    end
end
%**************************************************************************
%Plot of Lateral Constrain Inversion::
q = 2*l-1;
B = reshape(m_finforplot,q,n_d);

figure('visible', 'on')
v = [];
mag = [];
t_1 = (q+1)/2+1;
for i =1:n_d
    v = [v;i,0;i+1,0;i+1,-1*B(t_1,i);i,-1*B(t_1,i);];    
    t_2 = B(t_1,i);
    for j = 1:(q-1)/2-1
        v = [v;i,-1*t_2;i+1,-1*t_2;i+1,-1*(t_2+B(t_1+j,i));i,-1*(t_2+B(t_1+j,i));];
        t_2 = t_2+B(t_1+j,i);
    end
    v = [v;i,-1*t_2;i+1,-1*t_2;i+1,-250;i,-250;];
end

for i =1:n_d
    for j = 1 : (q+1)/2
        mag  = [mag; B(j,i)];
    end
end
V = v;
face = [1:4*n_d*((q+1)/2)];
face = reshape(face,4,(n_d*((q+1)/2)));
F = face';
C = mag;
patch('Faces',F,'Vertices',V,'FaceVertexCData',C,'FaceColor','flat');
colormap(parula);
colorbar;
str2fu = {'Res','Plot','Laterally_constrain'};
[path] = path_check(pwd,str2fu);
name = strcat('lC_',num2str(y+4),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');

disp('Plots and .dat file for every model is saved in Res Folder. ');
%************************************************************************
% Td_Res = reshape(Td_Res,301,n_d);
% 
% figure(12222)
% imagesc(Td_Res)
% map = interp1([0;1],[1 1 1;0 0.45 0.74],0:0.01:1);
% colormap(map)
% colorba
%************************************************************************
%RMSE plot

figure('visible', 'on');
plot(RMS,'-kx');
xlabel('Soundings');
ylabel('Root Mean Error');
title('RMSE Values Marq');

str2fu = {'Res','Plot','Laterally_constrain'};
[path] = path_check(pwd,str2fu);
name = strcat('lC_',num2str(y+2),'.jpg');
saveas(gcf, fullfile(path, name), 'jpg');
