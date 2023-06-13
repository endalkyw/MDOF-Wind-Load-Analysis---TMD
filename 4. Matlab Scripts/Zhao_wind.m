%Wind generation based on the Zhao et al. (2021) paper
function [x_j,u_bar]=Zhao_wind(exposure_category, mean_speed, heights, seed_value)
    global heights_;global friction_velocity_;global wind_velocities_;
    
    num_floors_=length(heights);
    w_cutoff=5*2*pi; %5Hz cutoff frequency
    if num_floors_==1
        M=6000; %to have a 10min time history
        N=3000;
    else
        M=4096;
        N=2048;
    end
    delta_w=w_cutoff/N;
    delta_t=(2*pi)/(M*delta_w);
    


    heights_=zeros(num_floors_,1);
    heights_(1)=heights(1);
    for i=2:length(heights)
        heights_(i)=heights_(i-1)+heights(i);
    end
    %Calculate the friction velocity and mean wind speeds based on the ASCE
    %exposure condition
    [friction_velocity_,wind_velocities_]=exposure_category_velocity(exposure_category,...
        heights_, 0.4,mean_speed);
    u_bar=wind_velocities_; %return the mean wind speeds at all floors

    %step (1) compute the Cholesky decomposition at each single indexing frequency
    %and generate independent random phase angles
    elements=(num_floors_^2+num_floors_)/2;
    H_vals=zeros(N,elements);
    for L=0:N-1
        S_L=cross_spectral_density((L*delta_w)+(delta_w/2));
        H_L=chol(S_L,'lower');
        count=1;
        for j=1:num_floors_
            for k=1:j
                H_vals(L+1,count)=H_L(j,k);
                count=count+1;
            end
        end
    end

    if strcmp(seed_value,'N')==1
        seed='shuffle'; %random number generation based on current time; gives different sequence at each call 
    else
        seed=str2double(seed_value);
    end
    
    %generating the independent uniformly distributed random phase angles
    rng(seed,'mt19937ar'); %use the default (Mersenne Twister) generator with the provided seed value
    phi_kL=2*pi*rand(N,num_floors_);
    x_j=zeros(M*num_floors_,num_floors_);
    j=1;
    count=1;
    %follow steps (2) through (5)
    while(j<num_floors_)
        %step (2) generate the frequency sequences F_jk on equation (29) and use
        %FFT to obtain the time sequences G_jk
        F_jk=zeros(M,j);
        for k=1:j
            F_jk(1:N,k)=2*sqrt(delta_w)*H_vals(:,count).*exp(complex(0,atan(imag(H_vals(:,count))./real(H_vals(:,count))))).*exp(complex(0,phi_kL(:,k)));            count=count+1;
        end
        G_jk=M*ifft(F_jk);
        %step (3) use equation (40) to calculate h_jq; use FFT to obtain
        %X_tilde_jq
        h_jq=zeros(M,num_floors_);
        q=(0:1:M-1)';
        for k=1:j
            h_jq(:,k+1)=G_jk(:,k).*exp(complex(0,k*(delta_w/num_floors_)*delta_t*q));%k+1 because the first column is reserved for k=0
        end
        %MATLAB ifft works column-wise so transpose h to do the ifft in the
        %num_floors_ dimension
        h_jq=h_jq';
        X_tilde_jq=num_floors_*ifft(h_jq);
        X_tilde_jq=X_tilde_jq';
        %step (4) generate the time history at the jth simulation point
        X_j=zeros(num_floors_*M,1);
        for i=1:num_floors_
            X_j((i-1)*M+1:i*M)=X_tilde_jq(:,i);
        end
        p=(0:1:(M*num_floors_)-1)';
        x_j(:,j)=real(X_j.*exp(complex(0,0.5*delta_w*delta_t.*p)));
        j=j+1;
    end
    %generate the wind at the last floor based on steps (6) and (7)
    F_nk=zeros(M,num_floors_);
    count=num_floors_-1;
    %On the first iteration of the following loop,elements-count gives first 
    %the index for the first entry of the last row. For e.g., for a 50 story building, 
    %H_vals(:,elements-(num_floors_-1)) gives the H(50,1) values for all
    %frequencies
    for k=1:num_floors_
        F_nk(1:N,k)=2*sqrt(delta_w)*H_vals(:,elements-count).*exp(complex(0,atan(imag(H_vals(:,elements-count))./real(H_vals(:,elements-count))))).*exp(complex(0,phi_kL(:,k)));        count=count-1;
    end
    G_nk=M*ifft(F_nk);
    h_nq=zeros(M,num_floors_);
    q=(0:1:M-1)';
    for k=1:num_floors_-1
        h_nq(:,k+1)=G_nk(:,k).*exp(complex(0,k*(delta_w/num_floors_)*delta_t.*q));
    end
    h_nq=h_nq'; %since ifft of a matrix is a column operation
    p_tilde=0:1:num_floors_-1; 
    X_tilde_nq=num_floors_*(ifft(h_nq))'+(G_nk(:,num_floors_).*exp(complex(0,delta_w*delta_t*q)).*exp(complex(0,delta_w*delta_t*p_tilde)));
    X_n=zeros(num_floors_*M,1);
    for i=1:num_floors_
        X_n((i-1)*M+1:i*M)=X_tilde_nq(:,i);
    end
    p=(0:1:(M*num_floors_)-1)';
    x_j(:,num_floors_)=real(X_n.*exp(complex(0,0.5*delta_w*delta_t.*p)));

end

function ret = cross_spectral_density(omega) 
    global friction_velocity_; global wind_velocities_; global heights_;  
    coherence_coeff = 10;%coefficient for coherence function
    cross_spectral_density = zeros(length(heights_),length(heights_));
    [row,cols] = size(cross_spectral_density);
    for i = 1 : row
        %remove(i) from friction velocity
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
