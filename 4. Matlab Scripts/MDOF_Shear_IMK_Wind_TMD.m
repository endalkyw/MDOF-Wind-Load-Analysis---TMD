function MDOF_Shear_IMK_Wind_TMD(dir)
    [G,B,A,E,TI] = getInputs("GI","BI","AI","EI","TI");
    
% -- 1. Import General Properties -----------------------------------------
    g=str2double(G(2));      % gravity
    StrengthLimitCheck=1;
    tol = 1e-4;              % str2double(G(3));

% -- 2. Import Analysis Properties ---------------------------------------------------
    uniformDampingRatio = str2double(A(3));
    pdeltaSwitch = str2double(A(2));  

    
% -- 3. Import Building Properties ----------------------------------------
    nStory=size(B,1);
    wi  = B(:,1);
    h   = B(:,2);
    k   = B(:,3);
    Fy  = B(:,4);
    a_s = B(:,5);
    dcdy = zeros(nStory,1);
    a_c= zeros(nStory,1);
    do = 0;
    Vo = 0;
    Xi = ones(nStory,1)*uniformDampingRatio;

    
% -- 4. Get wind forces ---------------------------------------------------
    Fwind = csvread('../3. Wind TH/WindLoad.csv');
    if(E(1)=="4")
        dt = 0.7854;
    else
        dt = 0.1;
    end
    
% -- 5. Get TMD properties ---------------------------------------------------
    % TI
  
% -- 4. Analyze -----------------------------------------------------------
    [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, TI, w] = MDOF_Wind_Analysis_TMD3(h, wi, pdeltaSwitch, k, Xi, Fwind, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, StrengthLimitCheck, TI);

%     out=[mean(u(4,2000:end)); max(abs(u(4,2000:end)-mean(u(4,2000:end)))); std(u(4,2000:end))]
%     time_=time';
%     u_=u(4,:)';
%     results_=[[0;0;0] out; time_ u_];
    
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
    
% -- 6. Frequency analysis of top floor displacements ---------------------------------------
%     dtmin=min(diff(time));
%     time_even=0:dtmin:time(end);
%     Xfft=interp1(time,(u(4,:)-mean(u(4,:))),time_even);
%         %Xfft=interp1(time,Fwind(5,:),time_even);
%         %Xfft=interp1(time,TO.RTMD(5,:),time_even);
%         %Xfft=interp1(time,(TO.uTMD-mean(TO.uTMD)),time_even);
%     Lfft=length(Xfft);
%     Ffft=1/dtmin;
%     Yfft = fft(Xfft);
%     P2fft = abs(Yfft/Lfft);
%     P1fft = P2fft(1:Lfft/2+1);
%     P1fft(2:end-1) = 2*P1fft(2:end-1);
%     f_fft = 2*pi*Ffft*(0:(Lfft/2))/Lfft;
%     semilogx(f_fft,P1fft)   
end