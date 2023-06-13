%Wind generation based on the Wittig & Sinha (1975) paper
function [wind,u_bar]=WS_wind(exposure_category, mean_speed, heights,seed_value, total_time,dt)
    global num_floors_; global seed_value_; global freq_cutoff_; global num_times_; 
    global num_freqs_; global frequencies_; global heights_; 
    global friction_velocity_; global wind_velocities_;
    
    num_floors_=length(heights);
    if strcmp(seed_value,'N')==1
        seed_value_=Inf;
    else
        seed_value_=str2double(seed_value);
    end
    freq_cutoff_=5.0; %in Hz
    time_step_=1.0/(2.0*freq_cutoff_);
    tmp1=total_time/time_step_;
    if mod(ceil(tmp1),2)==0
        num_times_=ceil(tmp1);
    else
        num_times_=ceil(tmp1)+1;
    end
    num_freqs_=num_times_/2;
    frequencies_=zeros(1,num_freqs_);
    for i=1: length(frequencies_)
        frequencies_(1,i)=i*(freq_cutoff_/num_freqs_);
    end
    %Calculate heights of each floor
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
    wind=generate();
end 

function wind_vels= generate()
    global heights_; global num_times_; 
    %Initialize wind velocity vectors
    wind_vels=zeros(num_times_,length(heights_));
    %complex_random_vals=zeros(num_freqs_, length(heights_));
    %Loop over heights to find time histories
    complex_random_vals=complex_random_numbers();
    for k=1: length(heights_)
        wind_vels(:,k)=gen_location_hist(complex_random_vals,k);
    end
end

function ret = cross_spectral_density(frequency)
    global heights_; global friction_velocity_; global wind_velocities_;
    %Coefficient for coherence function
    coherence_coeff = 10;

    cross_spectral_density = zeros(length(heights_),length(heights_));
    [row,cols] = size(cross_spectral_density);
    for i = 1 : row
        cross_spectral_density(i,i)=200.0*((friction_velocity_)^2)*heights_(i)...
            /(wind_velocities_(i)*(1.0+50.0*frequency*heights_(i)/wind_velocities_(i))^(5.0/3.0));
    end

    for i = 1 : row
        for j = i+1:cols
            cross_spectral_density(i,j) = sqrt(cross_spectral_density(i,i)...
                *cross_spectral_density(j,j))*exp(-coherence_coeff*frequency*abs(heights_(i)-heights_(j))/(0.5*(wind_velocities_(i) + wind_velocities_(j))));
        end
    end

    diag_vector = diag(cross_spectral_density);
    diag_mat = diag(diag_vector);
    ret = cross_spectral_density' + cross_spectral_density - diag_mat;
end

function ret = gen_location_hist(random_numbers,column_index)
    global num_freqs_;
    complex_full_range = zeros(1,2*num_freqs_);
    complex_full_range(2:1+num_freqs_)=random_numbers(1:num_freqs_,column_index);
    %two different ways to find the rest of complex_full_range values
    %complex_full_range(num_freqs_+2:num_freqs_+num_freqs_)=conj(flip(random_numbers(1:num_freqs_-1,column_index)));
    %complex_full_range(num_freqs_+1) = abs(random_numbers(num_freqs_,column_index));
    complex_full_range(num_freqs_+2: 2*num_freqs_) =conj(flip(complex_full_range(2:num_freqs_)));
    complex_full_range(num_freqs_+1)= abs(complex_full_range(num_freqs_+1));
    % Calculate wind speed using real portion of inverse Fast Fourier Transform
    % full range of random numbers
    node_time_history = ifft(complex_full_range);
    ret = node_time_history;
end

function complex_random=complex_random_numbers ()
    global seed_value_; global heights_; global num_freqs_; global frequencies_;
    global freq_cutoff_; global cross_spec_density_matrix;

    if seed_value_~=Inf
        seed=seed_value_+10;
    else
        %uses history seed (random number generation based on time)
        seed='shuffle';
    end
    rng(seed, 'mt19937ar');
    white_noise= zeros(length(heights_), num_freqs_);
    [row,col]=size(white_noise);
    for i=1:row
        for j=1:col
            white_noise(i,j)=complex(randn*sqrt(0.5),randn*imag(sqrt(complex(-0.5))));%or normrnd(0,sqrt(0.5))
        end 
    end
  %Iterator over all frequencies and generate complex random numbers
  %for discrete time series simulation
  %cross_spec_density_matrix=zeros(length(heights_), length(heights_));
  complex_random=zeros(num_freqs_, length(heights_));
  for i=1: length(frequencies_)
      %Calculate cross-spectral density matrix for current frequency
      cross_spec_density_matrix=cross_spectral_density(frequencies_(i));
      %Find lower Cholesky factorization of cross-spectral density
      lower_cholesky=chol(cross_spec_density_matrix,'lower');

      %This is Equation 5(a) from Wittig & Sinha (1975)
      complex_random(i,:)=num_freqs_* sqrt(2.0*freq_cutoff_/num_freqs_)*lower_cholesky*...
          white_noise(:,i); 
  end
end 
