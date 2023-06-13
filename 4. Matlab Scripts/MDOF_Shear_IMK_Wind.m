function MDOF_Shear_IMK_Wind(dir)
    [G,B,A,E] = getInputs("GI","BI","AI","EI"); clc;
     
% -- 1. Import General Properties -----------------------------------------
    pdeltaSwitch = str2double(A(2));  
    uniformDampingRatio = str2double(A(3));
    g=str2double(G(2));
    StrengthLimitCheck=1;
    tol = 10e-4;
    
% -- 2. Import Building Properties ----------------------------------------
    nStory=size(B,1);
    wi=B(:,1);                      % weight vectors
    h =B(:,2);                      % height ~
    k =B(:,3);                      % Stiffness ~ 
    Fy=B(:,4);                      % yield strength ~
    a_s=B(:,5);                     % hardning ratio ~
    dcdy = zeros(nStory,1);      % Vector of Duct. Capacity (dc/dy) 5
    a_c =  zeros(nStory,1);       % Post-capping Stiffness Ratio    -0.05
    do = 0;                         % initial displacement
    Vo = 0;                         % initial velocity
    Xi = uniformDampingRatio*ones(nStory,1);
    
% -- 3. Get wind forces ---------------------------------------------------
    Fwind = csvread('../3. Wind TH/WindLoad.csv');
    if(E(1)=="4")
        dt = 0.7854; % timestep for the chen and letchford downburst wind
    else
        dt = 0.1;    % timestep for the remaining wind models
    end
  
% -- 4. Analyze -----------------------------------------------------------
    [u, v, a_t, fs_st, fs, T, phi, M, K, C, time] = MDOF_Wind_Analysis(h, wi, pdeltaSwitch, k, Xi, Fwind, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, StrengthLimitCheck);

% -- 5. Save results ---------------------------------------------------  
    a_t = a_t*g;
    
    max_absolute_disp=zeros(nStory,1);
    rms_displacement=zeros(nStory,1);
    
    max_absolute_acce=zeros(nStory,1);
    rms_acceleration =zeros(nStory,1);
    
    max_absolute_vel=zeros(nStory,1);
    rms_vel =zeros(nStory,1);
    
    max_absolute_sf=zeros(nStory,1);
    rms_sf =zeros(nStory,1);
    
    for i=1:nStory
        max_absolute_disp(i,1) = max(abs(u(i,:)));
        rms_displacement(i,1)  = rms(u(i,1:end-1));
        
        max_absolute_acce(i,1) = max(abs(a_t(i,:)));
        rms_acceleration(i,1)  = rms(a_t(i,1:end-1));
        
        max_absolute_vel(i,1)=max(abs(v(i,:)));
        rms_vel(i,1) =rms(v(i,1:end-1));
        
        max_absolute_sf(i,1)=max(abs(fs_st(i,:)));
        rms_sf(i,1) =rms(fs_st(i,1:end-1));
    end
    x = getInputs("OC"); clc;
    
    writematrix(time,"..\2. Outputs"+dir+"\M1 TimeVector_TH.csv");
    
    if x(1)==1
        writematrix(u,"..\2. Outputs"+dir+"\M2 RelativeDisplacement_TH.csv");
    end
    if x(2)==1
        writematrix(fs_st,"..\2. Outputs"+dir+"\M3 ShearForce_TH.csv");
    end
    if x(3)==1
        writematrix(v,"..\2. Outputs"+dir+"\M4 Velocity_TH.csv");       
    end
    if x(4)==1
        writematrix(a_t,"..\2. Outputs"+dir+"\M5 Acceleration_TH.csv");
    end
    if x(5)==1
        writematrix(M,"..\2. Outputs"+dir+"\M6 MassMatrix.csv");
    end
    if x(6)==1
        writematrix(K,"..\2. Outputs"+dir+"\M7 StiffnessMatrix.csv");
    end
    if x(7)==1
        writematrix(phi,"..\2. Outputs"+dir+"\M8 ModeShape.csv");
    end
    if x(8)==1
        writematrix(C,"..\2. Outputs"+dir+"\M9 DampingMatrix.csv");
    end
    if x(9)==1
        writematrix(T,"..\2. Outputs"+dir+"\M10 FundamentalPeriod.csv");
    end
   
    if(x(10)==1)
        writematrix(max_absolute_disp,"..\2. Outputs"+dir+"\M11 MaxAbsDisplacement.csv");
    end
    if(x(11)==1)
        writematrix(rms_displacement,"..\2. Outputs"+dir+"\M12 rmsDisplacement.csv");
    end
    if(x(12)==1)
        writematrix(max_absolute_acce,"..\2. Outputs"+dir+"\M13 MaxAbsAcceleration.csv");
    end
    if(x(13)==1)
        writematrix(rms_acceleration,"..\2. Outputs"+dir+"\M14 rmsAcceleration.csv");
    end
    if(x(14)==1)
        writematrix(max_absolute_vel,"..\2. Outputs"+dir+"\M15 MaxAbsVelocity.csv");
    end
    if(x(15)==1)
        writematrix(rms_vel,"..\2. Outputs"+dir+"\M16 rmsVelocity.csv");
    end
    if(x(16)==1)  
        writematrix(max_absolute_sf,"..\2. Outputs"+dir+"\M17 MaxAbsShearForce.csv");
    end
    if(x(17)==1)
        writematrix(rms_sf,"..\2. Outputs"+dir+"\M18 rmsShearForce.csv");
    end
    
end
