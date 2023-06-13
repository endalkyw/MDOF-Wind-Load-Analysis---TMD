function generate_wind_load()
    [B,E] = getInputs("BI","EI"); clc;

    density      = 1.225;                      %in kg/m3
    drag_coeff   = str2double(E(2));
    category     = E(3);
    meanspeed    = str2double(E(4));
    width        = str2double(E(5));
    seed         = E(6);
    heights      = B(:,2); 
    tot_time     = str2double(E(7));
    if(E(1)=="4")
        time_step = 0.7854; % timestep for the chen and letchford downburst wind
    else
        time_step = 0.1;    % timestep for the remaining wind models
    end
    
    if(E(1,1)=='1')       % Witting Sinha
       [wind,u_bar]=WS_wind(category,meanspeed,heights,seed,tot_time,time_step);
       wind=wind';
       
    elseif(E(1,1)=='2')   % Deodatis 
       [wind,u_bar]=Deodatis_wind(category,meanspeed,heights,seed,tot_time,time_step);
       wind=wind';
    
    elseif(E(1,1)=='3')   % Zhao et al.
       [wind,u_bar]=Zhao_wind(category,meanspeed,heights,seed);
       wind=wind';
    elseif(E(1,1)=='4')
         [wind,~,~,~,t,~,~,~] = ChenLetchfordWind(heights,str2double(E(8)),str2double(E(9)),str2double(E(10)));
          u_bar = zeros(size(B,1),1);
          wind = wind';
%           win_len=ceil(600/0.785398163397448); % 0.78 is the time step of the wind
    end
    
    [num_floors, num_times]=size(wind);        
    wind_forces=zeros(num_floors,num_times);
    windTimeHistory = zeros(num_floors,num_times);
    windArea = zeros(num_floors,1);
    %%
    for i=1:num_floors
        if i~=num_floors
            windArea(i)=width*0.5*(heights(i)+heights(i+1));
        else
            windArea(i)=width*0.5*heights(i);
        end
    end
    
    for i=1:num_floors
        for j=1:num_times
            windTimeHistory(i,j) = wind(i,j)+u_bar(i,1);
            wind_forces(i,j)=0.5*density*drag_coeff*(windTimeHistory(i,j))^2*windArea(i); % made a formula correction here
        end
    end
    writematrix(wind_forces(:,1:tot_time*(1/time_step)),'..\3. Wind TH\WindLoad.csv')
    writematrix(windTimeHistory(:,1:tot_time*(1/time_step)),'..\3. Wind TH\WindSpeed.csv');
    
    for i = (1:num_floors)
        writematrix(wind_forces(i,tot_time*(1/time_step)),"..\3. Wind TH\WL_"+i+".csv")
    end
end