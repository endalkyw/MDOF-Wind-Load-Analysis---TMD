function [fric_vel, velocity_prof]=exposure_category_velocity(exposure_category, ...
    heights, karman_const, mean_speed)
    roughness_ht = 0.0;power_exponent = 0.0;power_factor = 0.0; 
    %power factor is used because the mean wind speed input is assumed to be for Category C.
    %the power factors for categories B and D were calculated such that the
    %wind speeds at Zg are equal. The values of Zg for the different categories
    %were taken from Table 26.11-1 from ASCE 7 standard (2017)
    
    %Set friction velocity and mean wind speeds based on exposure category
    %Uses both the logarithmic law and the power law
    
    %Category A is removed from recent editions of the ASCE standard
%     if strcmp(exposure_category,'A')==1
%         roughness_ht = 2.0;
%         power_factor = 0.48;
%         power_exponent = 1/3;
    if strcmp(exposure_category,'B')==1
        roughness_ht = 0.3;
        power_factor = 0.68;
        power_exponent = 1/4;
    elseif strcmp(exposure_category,'C')==1
        roughness_ht = 0.02;
        power_factor = 1;
        power_exponent = 1/6.5;
    elseif strcmp(exposure_category,'D')==1
        roughness_ht = 0.005;
        power_factor = 1.18;
        power_exponent = 1/9;
    else
        disp('category is not valid, please check input value');     
    end
    velocity_prof=zeros(length(heights),1);
    for i=1: length(heights)
        velocity_prof(i,1)=power_factor*mean_speed*power((heights(i,1)/10.0),power_exponent);
    end
    fric_vel=power_factor*mean_speed*karman_const/log(10.0/roughness_ht);
    turbulence_intensity=1/log(10/roughness_ht);%the turbulence intensity at 10m
    fileID=fopen('../2. Outputs/WindReport.txt','w');
    fprintf(fileID,'Mean wind speed at 10m (Cat. C)\t  %5.2fm/s\n',mean_speed);
    fprintf(fileID,'Turbulence intensity at 10m\t  %3.1f%%\n',turbulence_intensity*100);
    fprintf(fileID,'Friction velocity\t\t %5.2fm/s\n',fric_vel);
    fprintf(fileID,'Power law exponent\t\t  %4.3f\n',power_exponent);
    fprintf(fileID,'Roughness length\t\t  %4.3fm\n',roughness_ht);
    fclose(fileID);
end

