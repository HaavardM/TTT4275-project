%clear all

A = 1;
ph = pi/8;
f0 = 10^5;
omega_0 = 2*pi*f0;
Fs = 10^6;
T = 1/Fs;

N = 513;
M = 2^10;
SNR_db = 10;
noise_std = sqrt(A^2/10^(SNR_db/10));

steps = 10^5;
omega_0_est = zeros(steps, 1);
angle_est = zeros(steps,1);

signal = A*exp(1i*(omega_0*[1:N]*T + ph));
for j=1:steps
    noise = normrnd(0, noise_std, 1, N) + 1i*normrnd(0, noise_std, 1, N);
    x = signal + noise;

    x_fft = fft(x,M)/M;


    [peak, pos] = max(x_fft);

    omega_0_est(j) = 2*pi*(pos)/(M*T);
    
    % TODO: The signs of these two components must be corrected
    angle_est(j) = angle(exp(1i*omega_0_est(j)*best_n0(N)*T)*-peak);

    if mod(j, steps/100) == 0 && false
        fprintf('%i%%\n',100* j/steps)
    end
end
clc;
%% Analyze

omega_err = (omega_0 - omega_0_est);
omega_0_var = var(omega_err);
phase_err = (ph - angle_est);
phase_var = var(phase_err);


%% Fminify
% Minimize error of our estimate
func = @(w, p) sum(abs(x - A*exp(1i*(w*[1:N]*T + p))));
% Start with the estimate from the ffts
[vals, ~, exitflag, output] = fminsearch(@(input) func(input(1), input(2)), [mean(omega_0_est), mean(angle_est)]);

omega_mle = vals(1);
phase_mle = vals(2);

%% Printing
fprintf('---- Input Data ----\n');
fprintf('A:\t%f\tomega:\t%fHz\tPhase offset:\t%fpi\n', A, omega_0, ph/pi);
fprintf('Fs:\t%fHz\tSNR:\t%idB\n', Fs, SNR_db);
fprintf('N: %i\t FFT-Size: %i\t#Simulations: %i\n\n', N, M, steps);

fprintf('---- Estimated Values (FFT-estimate) ----')
fprintf('Estimated omega_0: %.0f (%.3f%% off)\n', mean(omega_0_est), 100*(mean(omega_0_est) - omega_0)/omega_0);
fprintf('Error variance: %f | CRLB: %f\n\n', omega_0_var, best_freq(noise_std, A, T, N))
fprintf('Estimated angle: %f pi (%.3f%% off)\n', mean(angle_est)/pi,100*(mean(angle_est) - ph)/ph)
fprintf('Error variance: %f | CRLB: %f\n\n',  phase_var, best_phase(noise_std, A, N))

fprintf('---- After Minimization ----\n');
fprintf('Estimated omega_0: %.0f (%.3f%% off)\n', omega_mle, 100*(omega_mle - omega_0)/omega_0);
fprintf('Estimated angle: %f pi (%.3f%% off)\n', phase_mle/pi,100*(phase_mle - ph)/ph)


%% Helper functions
function snr = SNR(A, sigma)
    snr = A^2/(2*sigma^2);
end

function omega = best_freq(sigma, A, T, N) 
omega = 12*sigma^2 /(A^2 * T^2 * N * (N^2 - 1));
end
function phase = best_phase(sigma, A, N) 
       P = make_P(N);
       Q = make_Q(N);
       n0 = best_n0(N);
       phase = var_phase(sigma, n0, N, P, Q, A);
end
function phase = var_phase(sigma, n0, N, P, Q, A)
    phase = 12*sigma^2*(n0^2*N + 2*n0*P + Q)/(A^2*N^2*(N^2 - 1));
end
function P = make_P(N)
    P = N*(N-1)/2;
end
function Q = make_Q(N)
    Q = N*(N-1)*(2*N-1)/6;
end
function n0 = best_n0(N)
    n0 = make_n0(make_P(N), N);
end
function n0 = make_n0(P, N)
    n0 = -P/N;
end
