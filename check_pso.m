clear all;
clc;
disp('determine number of layers:  from 3-5:  ');
n = input('Enter number of Layers: ');

fig_input = 1;


n_d = input('Enter total number of soundings to check:: ');
str_var = hem_;
%**************************************************************************
RMS = [];
depth_final=[];
Resistivity_final =[];
RMS_var = 0;
Min_RMS = 99999;
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
    d0 = [d0
          re
          Im];
       
    dp=forward_HEM(rho,d,h,r,f);
    m=PSOnew(m,d0,rho,h,r,f,d);
    m
    dp=forward_HEM(m(1:length(rho)),m(length(rho)+1:length(m)),h,r,f);
    t = sqrt(norm((d0-dp)./d0)^2/length(d0));


    % Plots  of Predicted real and Imaginary plot
    if(fig_input == 1)
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
        
        figure('Visible', 'on');
        
        plot(F_Resistivity,F_depth,'r');
        LineWidth = 1.4;
        xlabel('Resistivity');
        ylabel('Depth');
        xlim([min(F_Resistivity)-10 max(F_Resistivity)+10]);
        title({strcat('(RMSE =  ',num2str(123),')'),'Resistivity vs Depth'});

    end
end