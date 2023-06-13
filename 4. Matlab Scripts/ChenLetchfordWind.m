function [wind,kappa,U_bar,ft,t,Vc,u,V_tz] =ChenLetchfordWind(heights,vmax,mf,zmax)
    global V_tz; global t; global num_floors_; global heights_; global V_c; 
    
    num_floors_=length(heights);
    heights_=zeros(num_floors_,1);
    heights_(1)=heights(1);
    for i=2:length(heights)
        heights_(i)=heights_(i-1)+heights(i);
    end
    
    ft = timefunc;
    Vc = V_c;
    Vertical_vel= V_z(vmax,zmax);
    V_tz = ft(1)*Vertical_vel; %See the explanation after equation (9)
    U_bar = Vertical_vel * ft;
    U_bar=U_bar';
    kappa = Deodatis_wind;
    
%     a = 0.25 * U_bar;
    a = mf * U_bar;
    u = a .* kappa;
    wind = U_bar + u;
    
end

%Wood's vertical velocity profile, equation (5) Chen--Letchford (2004)
function [ver_vel]=V_z(vmax,zmax)
    global heights_; 
    V_max = vmax;
    delta = 5*zmax;
%     delta = 400;
    ver_vel=1.55*(heights_/delta).^(1/6).*(1-erf(0.7*heights_/delta))*V_max;
end

function f_t=timefunc()  
    global M; global w_cutoff; global N;  global delta_w; global delta_t;
    global num_floors_; global t; global V_c
    w_cutoff=4; %in rad/s
    N=2048;
    M=4096;
    %(22) and (50b) from Deodatis (1996) as can be seen in (63 a-b) 
    delta_w=w_cutoff/N;
    delta_t=(2*pi)/(M*(w_cutoff/N));
    
    t = 0:delta_t:(num_floors_*M-1)*delta_t;
    V_t = 12;
    e = 150;
    d0 = 3000;
    r_mag = sqrt((d0-V_t*t).^2+ e^2);
    f = @Vel_r;
    V_r = arrayfun(f,r_mag);
    v_x = V_r.*(d0-V_t*t)./r_mag+V_t; %x is the direction of storm motion
    v_y = V_r * e./r_mag;
    V_c = sqrt(v_x.^2 + v_y.^2);
    f_t = V_c/max(V_c);  %equation (9), Chen--Letchford(2004)
end

%radial velocity distribution, see equation (6)
function vel=Vel_r(r)
    Vr_max = 47;
    r_max = 1000;  
    R = 700; 
    if r<r_max
        vel=Vr_max*(r/r_max);
    else
        vel=Vr_max * exp(-((r-r_max)/R)^2);
    end
end

%Wind generation based on the Deodatis (1996) paper
function [f_i]=Deodatis_wind()
    global phi_L; global num_floors_; global N;
    
    %generating the independent uniformly distributed random phase angles described
    %just before equation (24)
    rng('shuffle','mt19937ar'); %use the default (Mersenne Twister) generator with the provided seed value
    phi_L=2*pi*rand(N,num_floors_); %or phi_L=unifrnd(0,2*pi,N,num_floors_);
    %Now generate the wind
    f_i=calc_f();
end

function f_i=calc_f()
    global num_floors_; global M; global delta_w; global delta_t;
    f_i=zeros(num_floors_*M,num_floors_);
    h=gen_h();
    for p=1:num_floors_*M-1
        count=1;
        for j=1:num_floors_
            for k=1:j
                f_i(p+1,j)=f_i(p+1,j)+real(h(p+1,count)*exp(complex(0,k*(delta_w/num_floors_)*p*delta_t)));
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
    %We only need the loop evaluation for L<N because of the condition in 
    %equation (54)
    for L=0:N-1
        count=1;
        S_L=cross_spectral_density((L*delta_w)+(m*delta_w/num_floors_));
        H_L=chol(S_L,'lower');    
        for j=m:num_floors_
                Hjm_L=H_L(j,m);
                thetajm_L=angle(H_L(j,m)); % Equation (16)
                
                %Equation (47) 
                Bjm_L(L+1,count)=2*abs(Hjm_L)*sqrt(delta_w)*exp(-complex(0,thetajm_L))*exp(complex(0,phi_L(L+1,m)));
            count=count+1;
        end
    end
end


function ret = cross_spectral_density(omega) 
    global V_tz; global heights_;
    coherence_coeff = 10.0; %coefficient for coherence function
    cross_spectral_density = zeros(length(heights_),length(heights_));
    [row,cols] = size(cross_spectral_density);
    for i = 1 : row
        cross_spectral_density(i,i)=(1/(4*pi))*(200/6)*heights_(i)...
            /(V_tz(i)*(1.0+(50.0*omega*heights_(i))/(2*pi*V_tz(i)))^(5/3));
    end

    for i = 1 : row
        for j = i+1:cols
            %Davenport coherence
            cross_spectral_density(i,j) = sqrt(cross_spectral_density(i,i)...
                *cross_spectral_density(j,j))*exp(-coherence_coeff*(omega/(2*pi))*abs(heights_(i)-heights_(j))/(0.5*(V_tz(i) + V_tz(j))));
            
            %Fully correlated
%           cross_spectral_density(i,j) = sqrt(cross_spectral_density(i,i)...
%             *cross_spectral_density(j,j))*0.99;

            %Uncorrelated
%             cross_spectral_density(i,j) = sqrt(cross_spectral_density(i,i)...
%              *cross_spectral_density(j,j))*1e-6;
        end
    end

    diag_vector = diag(cross_spectral_density);
    diag_mat = diag(diag_vector);
    ret = cross_spectral_density' + cross_spectral_density - diag_mat;
end