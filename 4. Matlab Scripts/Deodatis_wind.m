%Wind generation based on the Deodatis (1996) paper
function [f_i,u_bar]=Deodatis_wind(exposure_category,mean_speed,heights,seed_value,tot_time,dt)
    global friction_velocity_; global wind_velocities_; global heights_;
    global num_floors_; global phi_L; 
    global M; global w_cutoff; global N;  global delta_w; global delta_t;
    
    num_floors_=length(heights);
    w_cutoff=(1/(dt*2))*2*pi; %in rad/s      fc = 50hz, fs = 100hz,  delta_T = 0.01
   
    M=2^(nextpow2((tot_time/(dt*num_floors_)))); % to have a 10min time history
    if (M<4096)
        M = 4096;
    end
        N=M/2;
    
    %(22) and (50b) from the paper as can be seen in (63 a-b) 
    delta_w=w_cutoff/N;
    delta_t=(2*pi)/(M*(w_cutoff/N));
    %T_0=num_floors_*((2*pi)/delta_w);
    
    heights_=zeros(num_floors_,1);
    heights_(1)=heights(1);
    for i=2:length(heights)
        heights_(i)=heights_(i-1)+heights(i);
    end
    %Calculate the friction velocity and mean wind speeds based on the ASCE
    %exposure condition
    [friction_velocity_,wind_velocities_]=exposure_category_velocity(exposure_category,...
        heights_, 0.4,mean_speed);
    u_bar=wind_velocities_;
    
    if strcmp(seed_value,'N')==1
        seed='shuffle'; %random number generation based on current time; gives different sequence at each call 
    else
        seed=str2double(seed_value);
    end
    
    %generating the independent uniformly distributed random phase angles described
    %just before equation (24)
    rng(seed,'mt19937ar'); %use the default (Mersenne Twister) generator with the provided seed value
    phi_L=2*pi*rand(N,num_floors_); %or phi_L=unifrnd(0,2*pi,N,num_floors_);
    %Now generate the wind
    f_i=calc_f();
end

function f_i=calc_f()
    global M; global delta_w; global delta_t; global num_floors_; 
    f_i=zeros(num_floors_*M,num_floors_);
    h=gen_h();
    for i=1:num_floors_*M
        count=1;
        for k=1:num_floors_
            for j=1:k
                f_i(i,k)=f_i(i,k)+real(h(i,count)*exp(complex(0,j*(delta_w/num_floors_)*(i-1)*delta_t)));
                count=count+1;
            end
        end
    end  
end

function h=gen_h()
    global M; global num_floors_;
    h_terms=(num_floors_^2+num_floors_)/2;
    h=zeros(num_floors_*M,h_terms);
    for m=1:num_floors_
        %when generating h_jm, we can generate h's with the same m value together. 
        %For e.g., h_11,h_21,h_31... can be generated together. But, we want to use h's with 
        %the same j value together to generate f. For e.g., we use h_31,h_32, and h_33 to 
        %generate f3. Here we want to generate all the h's and then store them in the order
        %they will be used by calc_f(). That's why we need the variables init and jump.
        g_m=gen_g(m);
        init=(m^2+m)/2;
        jump=m-1;
        for j=1:num_floors_-m+1
            h(1:M,init)=g_m(:,j);
            jump=jump+1;
            init=init+jump;    
        end
    end
    for k=1:num_floors_-1
        h(k*M+1:(k+1)*M,:)=h(1:M,:);%represents equation (45) from the paper where the h values are repeated
    end 
end
function g_jm=gen_g(m)
    global M; global B_jm; 
    B_jm=genB(m);
    g_jm=M*ifft(B_jm); %equation (51) from the paper. The summation is implemented using the IFFT.
end

function Bjm_L=genB(m)
    global M; global delta_w; global phi_L; global N;global num_floors_; 
    Bjm_L=zeros(M,num_floors_-m+1);
    for L=0:M-1
        count=1;
        %We only need the function evaluation for L<N because of the
        %condition in equation (54)
        if L<N
            S_L=cross_spectral_density((L*delta_w)+(m*delta_w/num_floors_));
            H_L=chol(S_L,'lower');
        end
        for j=m:num_floors_
            if L<N
                Hjm_L=H_L(j,m);
                thetajm_L=atan(imag(Hjm_L)/real(Hjm_L)); % Equation (16)
                %(47) From the paper
                Bjm_L(L+1,count)=2*abs(Hjm_L)*sqrt(delta_w)*exp(-complex(0,thetajm_L))*exp(complex(0,phi_L(L+1,m)));
            else
                Bjm_L(L+1,count)=0; %the condition in equation (54)
            end
            count=count+1;
        end
    end
end

function ret = cross_spectral_density(omega) 
    global friction_velocity_; global wind_velocities_; global heights_;
    coherence_coeff = 10.0;%coefficient for coherence function
    cross_spectral_density = zeros(length(heights_),length(heights_));
    [row,cols] = size(cross_spectral_density);
    for i = 1 : row
        cross_spectral_density(i,i)=(1/(4*pi))*200.0*((friction_velocity_)^2)*heights_(i)...
            /(wind_velocities_(i)*(1.0+(50.0*omega*heights_(i))/(2*pi*wind_velocities_(i)))^(5/3));
    end

    for i = 1 : row
        for j = i+1:cols
            cross_spectral_density(i,j) = sqrt(cross_spectral_density(i,i)...
                *cross_spectral_density(j,j))*exp(-coherence_coeff*(omega/(2*pi))*abs(heights_(i)-heights_(j))/(0.5*(wind_velocities_(i) + wind_velocities_(j))));
        end
    end

    diag_vector = diag(cross_spectral_density);
    diag_mat = diag(diag_vector);
    ret = cross_spectral_density' + cross_spectral_density - diag_mat;
end

