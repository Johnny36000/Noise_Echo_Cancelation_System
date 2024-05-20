[x, fs] = audioread("Shazli_Weird1.wav");      % Read audio file and sample rate

%%%%%%%%%%%%          SYSTEM EQ        %%%%%%%%%%%%%%%
%%%%%%%%%%%%%% y(n) = x(n) + ax(n-N) %%%%%%%%%%%%%%%%%

n = 1:size(x, 1);   % Independant variable (varies from 1 to the length of x[n])
N = 6e3;            % Delay by experimentation

a = 0.7;            % Amplitude of echo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SIMPLE ECHO %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo = a * x(n(N+1:end)-N);     % Echo signal alone
echo = [zeros(N, 1); echo];     % Echo signal preceeded with 0s for dimension matching

y = x + echo;                   % Add echo to the original system

y_normalised = y/max(abs(y));   % Normalisation for the use of audio write as it accepts only [-1, 1] (limitation of audiowrite() function)

audiowrite("ShazliEchoSys.wav", y_normalised, fs);   % Write audio file


%%% CONV %%%

h = 1 * (n==1) + a * (n==N);            % Impulse response of a simple echo generator system

yconv = conv(x, h, "full");             % Convolution step (x[n] ∗ h[n] = y[n]) using convolution function

yconv_normalised = yconv/max(abs(yconv)); % Normalization for yconv
audiowrite("ShazliEchoConv.wav", yconv_normalised, fs);  % Write audio file

%%% DECONV using FT %%% 

% Using the FT property F[n] ∗ G[n] = F(w)G(w)
% Use this to deconvolve H from D (delta) and get H's inverse

delta = 1 * (n==1);

D = fft(delta);             % Use the FFT (Discrete FT) algo to do FT
H = fft(h);
H2 = D ./ H;
h2 = ifft(H2);              % Use the IFFT (Inverse FFT) to go back to time domain

clean = conv(yconv, h2, "full");        % Apply convolution to the signal with echo with the inverse response to remove echo

clean_normalised = clean/max(abs(clean));
audiowrite("cleanS.wav", clean_normalised, fs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% COMPLEX ECHO %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

N = 14e3;               % Changed delay and amplitude of the echo to make it clearer
a = 0.5;

compEcho = zeros(n(end), 1);
for i = 1:size(n, 1) / N
    index_range = (i-1)*N + 1: i*N; 

    update_dim = max(index_range) + numel(a) * (i-1); 

    compEcho(1:update_dim) = compEcho(1:update_dim) + [zeros(i*N, 1); a.^i .* x(n(i*(index_range)) - i*N)];   % Generate echo ERROR: arrays are not the same size 
end

yComp = x + compEcho;       % Add echo to 
% the original signal

yComp_normalised = yComp/max(abs(yComp));
audiowrite("ShazliCECHOsys.wav", yComp_normalised, fs);      % Normalise and write the file

%%%%% CONV %%%%%

h11 = 1 * (n==1);                           % Generate impulse response for infinite echo
for i = 1:9                                 % With trial and error, 9 reflections turned out to be a sweet spot
    h11 = h11 + a^i .* (n == i*N);
end

yCompConv = conv(x, h11, "full");

yCompConv_normalised = yCompConv/max(abs(yCompConv));
audiowrite("ShazliCECHOconv.wav", yCompConv_normalised, fs);

%%%%% DECONV %%%%%
% Use the same FT property mentioned above

H11 = fft(h11);
H22 = D ./ H11;
h22 = ifft(H22); 

cleanComp = conv(yCompConv, h22, "full");
cleanComp_normalised = cleanComp/max(abs(cleanComp));
audiowrite("cleanComp.wav", cleanComp_normalised, fs);
tau= 1:2*n(end)+1;
figure(1);
plot(n,x);
title('Original sound');
grid;
figure(1);
plot(n,y);
title('Simple echo using system equation');


figure(2);
plot(tau,yconv);
title('Simple echo using convolution');
grid;

tau = 1:length(clean_normalised);

figure(3);
plot(tau,clean_normalised);
title('Simple echo cancelled');
grid;

figure(4);
plot(n,h);
title('Impulse response');
grid;


figure(5);
plot(n,h2);
title('impulse response of echo cancellation');
grid;


figure(6);
plot(n,yComp);
title('Infinite echo using system equation');
grid;

tau = (0:length(yCompConv)-1) / fs;  % Create a time axis based on the length of yCompConv and the sampling frequency

figure(7);
plot(tau, yCompConv);
title('Infinite echo using convolution');
xlabel('Time (s)');
ylabel('Amplitude');
grid;

tau = (0:length(cleanComp)-1) / fs;  % Create a time axis based on the length of cleanComp and the sampling frequency

figure(8);
plot(tau, cleanComp);
title('Infinite echo cancelled');
xlabel('Time (s)');
ylabel('Amplitude');
grid;