function [u, v, a_t, fs_st, fs, T, phi, M, K, C, time, TI, w] = ...
    MDOF_Wind_Analysis(h, wi, pdswitch, k, Xi, Fwind, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, StrengthLimitCheck, TI)


% Flipped sign of a


% [u, v, a_t, fs_st, fs, T, phi, M, K, C, time] = MDOF_ShearBiLinear_seismic(h, wi, Pdswitch, k, Xi, ug, dt, do, Vo, Fy, a_s, dcdy, a_c, tol, g, Name, StrengthLimitCheck)
% 
% Computes the WIND-DRIVEN response of a nonlinear Ibarra-Medina-Krawinkler 
% MDOF system, using Newmark's method, with Newton-Raphson iterations. 
% (Chopra 4th Edition, Table 16.3.3). The MDOF system is considered as a 
% shear building. It uses constant (average) acceleration method.
% 
% Inputs:
%   - h : Vector with story heights, in [length].
%   - wi: Vector with story weights, from 1 to the roof (end of the
%         vector), in [force].
%   - Pdswitch: If Pdelta = 0, no P-delta effects are considered. If Pdelta = 1 P-delta is considered.***
%   - k: Vector with story stiffnesses, from 1 to the roof (end of the
%        vector), in [force/length].
%   - Xi: Vector with damping ratio of each mode.
%   - Fwind: Matrix containing the time history of wind loads on the structure, in [N].
%           Rows are the stories, columns are the time step
%   - dt: Time step. Recommendation: Use dt = min[dt_crit,T/10,sampling of Fwind]
%           Where   dt_crit = inf ; for Average Acceleration method
%                   dt_crir = 0.551*dcmin(T) ; for Linear Acceleration method
%   - do: Vector with initial displacement, in [length].
%   - Vo: Vector with initial velocity, in [length/sec].
%   - Fy: Vector of Yield Force of each story i, in [force].
%   - a_s: Vector of hardening stiffness ratio of each story.
%   - dcdy: Vector of ductility capacity (dc/dy) for capping point of each story.
%   - a_c: Vector of post-capping stiffness ratio of each story.
%   - tol: Tolerance of the Newton-Raphson iterations. 
%           Recommendation: Use tol = 1E-3 to 1E-8.
%   - g: Standard gravity acceleration, in [length/sec^2].
%   - Name: Name of the ground motion [string].
%   - StrengthLimitCheck:   1 if strength limit is considered (recommended)
%                           0 if strength limit is NOT considered
%   - TI = TMD Inputs
% 
% Outputs:
%   - u: Matrix with Relative Displacement time history, in [length].
%   - v: Matrix with Relative Velocity time history, in [length/sec].
%   - a_t: Matrix with Absolute Acceleration time history, in [g].
%   - fs_st: Matrix with Story Restoring Force (Shear Force) time history, in [force].
%   - fs: Matrix with Story Forces time history, in [force].
%   - T: Vector of periods of the structure, in [sec].
%   - phi: Matrix with modal shapes of the structure.
%   - M, K, C: Mass, Stiffness and Damping matrices assambled, with
%              consistent units.
%   - time: Time vector, in [sec].
%
%   u, v, a_t and Fs, have N rows, with N being equal to the number of
%   stories (N = length(M) = length(K)), and the same number of columns
%   as Fwind, representing the time variation of the response.
%
%   All the units must be consistent between the input parameters. The
%   output parameters will be also consistent.
%  
%   Sesimic version by Pablo Heresi (08-08-2015)
%   Wind version by Eric Lalonde (06-12-2022)


Pi = determine_Pi(wi,pdswitch);
N = length(wi);     % Number of stories
MaxIter = 20;

%% Obtain T and phi

% Note: first index corresponds to 1st floor, and last index to roof.
% M and K matrices
M = diag(wi)/g;         % Mass Matrix
K = ComputeK1(k);        % Stiffness Matrix

% Eigenvalue analysis
[phi,w2] = eig(K,M);
w = sqrt(diag(w2));     % Undamped frequencies
[w,index] = sort(w);
T = 2*pi./w;            % Undamped periods

% Save initial frequencies without TMD
phi_init=phi;
w_init=w;

% Initial Damping Values
% Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
sphi = phi;
for i = 1:N
    sphi(:,i) = phi(:,index(i))/ max(abs(phi(:,index(i))));
end
phi = sphi;             % Normalized modal shapes

% C matrix
Mi = diag(phi'*M*phi);
if size(Xi) == size(Mi)
    Ci = 2*Mi.*w.*Xi;
else
    Ci = 2*Mi.*w.*Xi';
end
C = (phi')^(-1)*diag(Ci)*phi^(-1);

% % C Matrix using Rayleigh Damping
% f1=w(1)/(2*pi);
% E1=Xi(1);
% if N>1
%     f2=w(2)/(2*pi);
%     E2=Xi(2);
% else
%     f2=f1;
%     E2=E1;
% end
% 
% alpha=4*pi*f1*f2*(E1*f2-E2*f1)/(f2^2-f1^2);
% beta=(E2*f2-E1*f1)/(pi*(f2^2-f1^2));
% C=alpha*M+beta*K;



% If there are TMDs
if TI.nTMD>0

    % Get TMD properties
    % Get target frequency of each TMD
    TI.freqTMD=zeros(TI.nTMD,1);
    for iTMD=1:TI.nTMD
        if isnumeric(TI.targFreqTMD{iTMD}) % If a specific frequency is specified
            TI.freqTMD(iTMD)=TI.targFreqTMD{iTMD};
        elseif strcmp(TI.targFreqTMD{iTMD}(1),'M') % If "Mode 4" is specifified, e.g.
            modeNumTMD=str2double(erase(TI.targFreqTMD{iTMD},'Mode '));
            modeFreqTMD=w(modeNumTMD); % Get frequency of target mode
            freqFactorTMD=1-1.2*TI.massRatioTMD(iTMD); % Qj factor from Den Hartog
            TI.freqTMD(iTMD)=freqFactorTMD*modeFreqTMD;
        end
    end

    % Optimal damping ratios (Xi_j in Den Hartog)
    TI.dampRatioTMD=((TI.massRatioTMD.*(3-(0.5.*TI.massRatioTMD).^0.5)) /...
        (8.*(1+TI.massRatioTMD).*(1-0.5.*TI.massRatioTMD))).^0.5;

    % Mass
    TI.mTMD=TI.massRatioTMD.*wi(TI.locTMD)/g;

    % Stiffness
    TI.kTMD=TI.mTMD.*TI.freqTMD.^2;


    % Add TMD properties to structural matrices
    % Mass
    M_init=M;
    M(N+(1:TI.nTMD),N+(1:TI.nTMD))=zeros(TI.nTMD); % Add zero rows and columns to end
    for iTMD=1:TI.nTMD
        M(N+iTMD,N+iTMD)=TI.mTMD(iTMD); % Add TMD masses to diagonal of mass matrix
    end
    
    % Stiffness
    K_init=K;
    K(N+(1:TI.nTMD),N+(1:TI.nTMD))=zeros(TI.nTMD); % Add zero rows and columns to end
    for iTMD=1:TI.nTMD
        K(N+iTMD,N+iTMD)=TI.kTMD(iTMD);  % Add TMD stiffnesses to diagonal of stiffness matrix
        K(TI.locTMD(iTMD),TI.locTMD(iTMD))=K(TI.locTMD(iTMD),TI.locTMD(iTMD))+TI.kTMD(iTMD);  
            % Add TMD stiffnesses to diagonal of original structural matrix
        % Add negtive TMD stiffnesses to intersection of story DOF and TMD DOF
        K(TI.locTMD(iTMD),N+iTMD)=-TI.kTMD(iTMD);
        K(N+iTMD,TI.locTMD(iTMD))=-TI.kTMD(iTMD); 
    end
    
    % Mass for P-Delta
    TI.wi=wi;
    for iTMD=1:TI.nTMD
        TI.wi(TI.locTMD(iTMD))=TI.wi(TI.locTMD(iTMD))+g*TI.mTMD(iTMD);
    end
    Pi = determine_Pi(TI.wi,pdswitch);
    
    
    % Add empty cells to force time history
    Fwind_init=Fwind;
    Fwind(N+(1:TI.nTMD),:)=zeros(TI.nTMD,length(Fwind(1,:))); % Add zero rows to bottom
    
    % Add empty cells to additional properties
    % Pi
    Pi_init=Pi;
    Pi(N+(1:TI.nTMD))=zeros(TI.nTMD,1);
    % h
    h_init=h;
    h(N+(1:TI.nTMD))=ones(TI.nTMD,1);
    
    
    % Rerun eigenvector analysis
    % Eigenvalue analysis
    [phi,w2] = eig(K,M);
    w = sqrt(diag(w2));     % Undamped frequencies
    [w,index] = sort(w);
    T = 2*pi./w;            % Undamped periods
    
    
    
    % Damping
    Xi_init=Xi;
    Xi(N+(1:TI.nTMD))=zeros(TI.nTMD,1); % Add zero rows to end
    for iTMD=1:TI.nTMD
        Xi(N+iTMD)=TI.dampRatioTMD(iTMD); % Add TMD damping ratios to damping array
    end
    
    % Find where TMD damping needs to be arranged in global Xi matrix
    [~,freqSort]=sort([w_init;TI.freqTMD]);
    Xi=Xi(freqSort);
    for iTMD=1:TI.nTMD
       TI.freqSort(iTMD)=find(freqSort==(N+iTMD)); 
    end
    
    % Calc C of different TMDs
    TI.cTMD=2.*TI.mTMD.*TI.dampRatioTMD.*TI.freqTMD;
    
    % Add C value for TMD on diagonals of originally developed C matrix
    C_init=C;
    C(N+(1:TI.nTMD),N+(1:TI.nTMD))=zeros(TI.nTMD); % Add zero rows and columns to end
    for iTMD=1:TI.nTMD
        C(N+iTMD,N+iTMD)=TI.cTMD(iTMD);  % Add TMD damps to diagonal of damp matrix
        C(TI.locTMD(iTMD),TI.locTMD(iTMD))=C(TI.locTMD(iTMD),TI.locTMD(iTMD))+TI.cTMD(iTMD);  
            % Add TMD damps to diagonal of original structural matrix
        % Add negtive TMD damp to intersection of story DOF and TMD DOF
        C(TI.locTMD(iTMD),N+iTMD)=-TI.cTMD(iTMD);
        C(N+iTMD,TI.locTMD(iTMD))=-TI.cTMD(iTMD); 
    end
    
    
    
    % Adjust N from number of stories to total number of DOF
    N_init=N;
    N=N+TI.nTMD;
    TI.N_init=N_init;
end
    TI.N=N;  
    

% Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
sphi = phi;
for i = 1:N
    sphi(:,i) = phi(:,index(i))/ max(abs(phi(:,index(i))));
end
phi = sphi;             % Normalized modal shapes



%% Check stability of the method

Np = length(Fwind);     % Length of the record
time = 0:dt:(dt*(Np-1));

% If the time step is too large, display a warning.
if dt > T(1)/20
    Warning = "Warning: The time step used (dt) for the wind time history is \ngreater than T1/20. This is not recommended for representing \nthe time response correctly.\n";
    fprintf(Warning)
    disp(['A new dt = T/20 = ' num2str(T(1)/20) ' sec is used.']);
   
    dt_ = T(1)/20;
    
    time_ = 0:dt_:600-dt_;
    
    Fwind_ = zeros(N,length(time_));

    for i=1:N        
          Fwind_(i,:) = interp1(time,Fwind(i,:),time_);
    end
    time = time_;
    dt = dt_;

    Fwind = Fwind_;
    
    Np = length(Fwind);
    
end



%% Initial Calculations

r = ones(N,1);          % Note that this assumes horizontal excitation
P = Fwind;              % Wind load

dy = Fy./k;             % Yielding displacement of each story
dc = dcdy.*dy;          % Capping displacement of each story
Fmax = (1-a_s).*Fy+a_s.*k.*dc;      % Positive Strength Limit
Fmin = -Fmax;                       % Negative Strength Limit
LimMax = zeros(N,1);
LimMin = zeros(N,1);

% Initialize vectors
fs_st = zeros(N,Np);        % Story restoring force
fs = fs_st;                 % Total floor restoring force
fs(:,1) = [-diff(fs_st(:,1));fs_st(end,1)];
Kt = K;                     % Initial tangent stiffness
kt = k;

u = zeros(N,Np);            % Relative displacement time history
v = u;                      % Relative velocity time history
a = u;                      % Relative acceleration time history
u(:,1) = do;                % Initial Displacement
v(:,1) = Vo;                % Initial Velocity
a(:,1) = M\(P(:,1)-C*v(:,1)-fs(:,1));  % Initial Relative Acceleration

% Constants
a1 = 4/dt^2*M + 2/dt*C;
a2 = 4/dt*M + C;
dt2 = dt/10;
a1_2 = 4/dt2^2*M + 2/dt2*C;
a2_2 = 4/dt2*M + C;

R = zeros(N,Np);            % Unbalanced force history
It = zeros(1,Np);           % Number of iterations
It(1) = 1;

%% Calculation for each time step
kt_prev = k;
i = 2;
kk = 0;
checkpoint=1000;

while i < size(u,2)
    
    if i > checkpoint
        ['Step ',num2str(i),' of ',num2str(length(time))...
            ,'. ',num2str(length(time)-i+1),' remaining'] %#ok<NOPRT>
        checkpoint=checkpoint+1000;
    end

    % Initial estimate of structure displacement and velocity
    u(:,i+1) = u(:,i);
    v(:,i+1) = 2/dt*(u(:,i+1)-u(:,i)) - v(:,i);

    % Initial estimate of internal and external forces
    fs_st(:,i+1) = fs_st(:,i);
    fs(:,i+1) = fs(:,i);

    p_ = P(:,i+1) + a1*u(:,i) + a2*v(:,i) + M*a(:,i);
    
    
    % Initial estimate of TMD displacement and velocity:
    currentDT=time(i+1)-time(i);
    
    % Determine accelerations
    a(:,i+1) = M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1));
    
    % Newton-Raphson iterations
    j = 0;
    R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1);
    
    while sum(abs(R(:,i+1)) > tol)   &&  j < MaxIter+1
        
        % Estimate structure displacement and velocity
        Kt_ = Kt + a1;
        du = Kt_\R(:,i+1);
        u(:,i+1) = u(:,i+1) + du;
        v(:,i+1) = 2/dt*(u(:,i+1)-u(:,i)) - v(:,i);
        
        % Estimate TMD structure and displacement
        currentDT=time(i+1)-time(i);
        
        % State determination
        [fs_st(:,i+1),kt,Fmax,Fmin,LimMax,LimMin] =...
            StateDet(h,Pi,u(:,max(i-1,1)),u(:,i),u(:,i+1),fs_st(:,i),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck,TI);
        if TI.nTMD>0
            % Internal story forces
            fs(1:TI.N_init,i+1) = [-diff(fs_st(1:TI.N_init,i+1));fs_st(TI.N_init,i+1)];
            % TMD forces
            fs((TI.N_init+1):end,i+1) = fs_st((TI.N_init+1):end,i+1);
            % Subtract TMD forces from stories
            for iTMD=1:TI.nTMD
                fs(TI.locTMD(iTMD),i+1)=fs(TI.locTMD(iTMD),i+1)-fs(TI.N_init+iTMD,i+1);
            end
        else
            fs(:,i+1) = [-diff(fs_st(:,i+1));fs_st(end,i+1)];
        end
        
        Kt = ComputeK2(kt,TI);
        
        % Determine accelerations
        a(:,i+1) = M\(P(:,i+1)-C*v(:,i+1)-fs(:,i+1));

        R(:,i+1) = p_ - fs(:,i+1) - a1*u(:,i+1);      % Unbalanced force
        j = j+1;                                % Increase # of iterations

    end
    if [kt == kt_prev  ;  j < MaxIter+1]
        kk = kk+1;
        It(i+1) = j;    % Total number of iterations
        kt_prev = kt;
        i = i+1;
    else    % Change in any stiffness or convergence not achieved

        windInsert = (Fwind(:,i+1)-Fwind(:,i))/10;
 
        time_int = (time(i)+dt2):dt2:(time(i+1)-dt2);
        
        windInsertMatrix=ones(N,9);
        for j=1:N
            if windInsert(j) == 0
                windInsertMatrix(j,:)=Fwind(j,i)*windInsertMatrix(j,:);
            else 
                windInsertMatrix(j,:)=Fwind(j,i)+windInsert(j)*windInsertMatrix(j,:);
            end
        end           
        
        time = [time(:,1:i) time_int time(:,i+1:end)];
        
        Fwind = [Fwind(:,1:i) windInsertMatrix Fwind(:,i+1:end)];

        P = Fwind;          % Equivalent external load

        u = [u(:,1:i) zeros(N,9) u(:,i+1:end)];

        fs_st = [fs_st(:,1:i) zeros(N,9) fs_st(:,i+1:end)];
        fs = [fs(:,1:i) zeros(N,9) fs(:,i+1:end)];
        v = [v(:,1:i) zeros(N,9) v(:,i+1:end)];
        a = [a(:,1:i) zeros(N,9) a(:,i+1:end)];
        R = [R(:,1:i) zeros(N,9) R(:,i+1:end)];
        It = [It(:,1:i) zeros(1,9) It(:,i+1:end)];
        
        
        for i2 = 1:10
        	u(:,i+i2) = u(:,i+i2-1);
            v(:,i+i2) = 2/dt2*(u(:,i+i2)-u(:,i+i2-1)) - v(:,i+i2-1);
            fs_st(:,i+i2) = fs_st(:,i+i2-1);
            fs(:,i+i2) = fs(:,i+i2-1);
            
          
            % Determine accelerations
            a(:,i+i2) = M\(P(:,i+i2)-C*v(:,i+i2)-fs(:,i+i2));
            
             
            p_ = P(:,i+i2) + a1_2*u(:,i+i2-1) + a2_2*v(:,i+i2-1) + M*a(:,i+i2-1);
                
            % Newton-Raphson iterations
            j = 0;
            R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2);
                
            while sum(abs(R(:,i+i2)) > tol)   &&  j < MaxIter+1
                Kt_ = Kt + a1_2;
                du = Kt_\R(:,i+i2);
                u(:,i+i2) = u(:,i+i2) + du;
                v(:,i+i2) = 2/dt2*(u(:,i+i2)-u(:,i+i2-1)) - v(:,i+i2-1);


                % State determination
                [fs_st(:,i+i2),kt,Fmax,Fmin,LimMax,LimMin] = StateDet...
                    (h,Pi,u(:,i+i2-2),u(:,i+i2-1),u(:,i+i2),fs_st(:,i+i2-1),k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck,TI);
                if TI.nTMD>0
                    % Internal story forces
                    fs(1:TI.N_init,i+i2) = [-diff(fs_st(1:TI.N_init,i+i2));fs_st(TI.N_init,i+i2)];
                    % TMD forces
                    fs((TI.N_init+1):end,i+i2) = fs_st((TI.N_init+1):end,i+i2);
                    % Subtract TMD forces from stories
                    for iTMD=1:TI.nTMD
                        fs(TI.locTMD(iTMD),i+i2)=fs(TI.locTMD(iTMD),i+i2)-fs(TI.N_init+iTMD,i+i2);
                    end
                else
                    fs(:,i+i2) = [-diff(fs_st(:,i+i2));fs_st(end,i+i2)];
                end

                Kt = ComputeK2(kt,TI);

                % Determine accelerations
                a(:,i+i2) = M\(P(:,i+i2)-C*v(:,i+i2)-fs(:,i+i2));

                R(:,i+i2) = p_ - fs(:,i+i2) - a1_2*u(:,i+i2);      % Unbalanced force
                j = j+1;                                % Increase # of iterations

            end
            It(i+i2) = j;    % Total number of iterations
        end
        i = i+10;
        kt_prev = kt;
    end

end
a_t = a/g; 	% Absolute acceleration, in [g]

end

function [K] = ComputeK1(k)

if length(k) > 1
    k_aux = k(2:end);
    k_aux(end+1,1) = 0;
    K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
else
    K = k;
end
end

function [K] = ComputeK2(k,TI)

if TI.N > 1
    if TI.nTMD>0
        k_aux = k(2:end);
        k_aux(end+1,1) = 0;
        
        if TI.N_init>1
            K(1:TI.N_init,1:TI.N_init) = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
        else
            K=k;
        end
        
        % TMD Stuff
        for iTMD=1:TI.nTMD
            K(TI.N_init+iTMD,TI.N_init+iTMD)=TI.kTMD(iTMD);  % Add TMD stiffnesses to diagonal of stiffness matrix
            K(TI.locTMD(iTMD),TI.locTMD(iTMD))=K(TI.locTMD(iTMD),TI.locTMD(iTMD))+TI.kTMD(iTMD);  
                % Add TMD stiffnesses to diagonal of original structural matrix
            % Add negtive TMD stiffnesses to intersection of story DOF and TMD DOF
            K(TI.locTMD(iTMD),TI.N_init+iTMD)=-TI.kTMD(iTMD);
            K(TI.N_init+iTMD,TI.locTMD(iTMD))=-TI.kTMD(iTMD); 
        end
    else
        k_aux = k(2:end);
        k_aux(end+1,1) = 0;
        K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);
    end
else
    K = k;
end
end
 
function [fs,kt,Fmax2,Fmin2,LimMax2,LimMin2] = StateDet(h,Pi,u0,u1,u2,fs1,k,Fy,a_s,dc,a_c,Fmax,Fmin,LimMax,LimMin,StrengthLimitCheck,TI)
if TI.nTMD>0
    du0 = [u0(1);diff(u0(1:TI.N_init));zeros(TI.nTMD,1)];
    du1 = [u1(1);diff(u1(1:TI.N_init));zeros(TI.nTMD,1)];
    du2 = [u2(1);diff(u2(1:TI.N_init));zeros(TI.nTMD,1)];
    for iTMD=1:TI.nTMD
       du0(TI.N_init+iTMD)=u0(TI.N_init+iTMD)-u0(TI.locTMD(iTMD));
       du1(TI.N_init+iTMD)=u1(TI.N_init+iTMD)-u1(TI.locTMD(iTMD));
       du2(TI.N_init+iTMD)=u2(TI.N_init+iTMD)-u2(TI.locTMD(iTMD));
    end
else
    du0 = [u0(1);diff(u0)];
    du1 = [u1(1);diff(u1)];
    du2 = [u2(1);diff(u2)];
end

Fmax2 = Fmax;
Fmin2 = Fmin;
LimMax2 = zeros(length(du1),1);
LimMin2 = zeros(length(du1),1);

fs = zeros(size(fs1));
fs1 = fs1 + Pi./h.*du1;
kt = k;
if TI.nTMD>0
    ilim=TI.N_init;
    for i=(TI.N_init+1):TI.N % Skip the below ductility checks
        fs(i)=fs1(i) + TI.kTMD(i-TI.N_init)*(du2(i)-du1(i));  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end
else
    ilim=TI.N;
end
for i = 1:ilim
    fs(i) = fs1(i) + k(i)*(du2(i)-du1(i));
    
    if du2(i) > dc(i)
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)-dc(i));
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);
        ktEnv = a_c(i)*k(i);
    elseif du2(i) > -dc(i)
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);
        NegLimEnv = (a_s(i)-1)*Fy(i)+a_s(i)*k(i)*du2(i);
        ktEnv = a_s(i)*k(i);
    else
        PosLimEnv = (1-a_s(i))*Fy(i)+a_s(i)*k(i)*du2(i);
        NegLimEnv = (a_s(i)-1)*Fy(i)-a_s(i)*k(i)*dc(i)+a_c(i)*k(i)*(du2(i)+dc(i));
        ktEnv = a_c(i)*k(i);
    end

    if fs(i) > min(Fmax(i),PosLimEnv)
        LimMax2(i) = 1;
        fs(i) = min(Fmax(i),PosLimEnv);
        if fs(i) == Fmax(i)
            kt(i) = 0;
        else
            kt(i) = ktEnv;
        end
    elseif fs(i) < max(Fmin(i),NegLimEnv)
        LimMin2(i) = 1;
        fs(i) = max(Fmin(i),NegLimEnv);
        if fs(i) == Fmin(i)
            kt(i) = 0;
        else
            kt(i) = ktEnv;
        end
    end
    
    if StrengthLimitCheck
        if du1(i) > dc(i) && du0(i) < du1(i) && du2(i) < du1(i) && LimMax(i)
            Fmax2(i) = fs1(i);
        end
        
        if du1(i) < -dc(i) && du0(i) > du1(i) && du2(i) > du1(i) && LimMin(i)
            Fmin2(i) = fs1(i);
        end
    end
    
    fs(i) = fs(i) - Pi(i)/h(i)*du2(i);
end                          
end

function [PI] = determine_Pi(w, ps)
  P = zeros(length(w),1);
  if(ps==1)
      P=w;
      for i = length(w)-1:-1:1
              P(i,1)=P(i,1)+P(i+1,1);
      end
  end
  PI = P;
end

