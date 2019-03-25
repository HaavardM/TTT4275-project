clear all

A = 1;
phase = pi/8;
f0 = 10^5;
omega_0 = 2*pi*f0;
Fs = 10^6;
T = 1/Fs;

N = 513;
M = 2^10;
SNR_db = 10;
noise_std = sqrt(A^2/10^(SNR_db/10));

steps = 10^4;
omega_0_est = zeros(1,steps);
angle_est = zeros(1, steps);

signal = A*exp(1i*(omega_0*[1:N]*T + phase));
for j=1:steps
    noise = normrnd(0, noise_std, 1, N) + 1i*normrnd(0, noise_std, 1, N);
    x = signal ;%+ noise;

    x_fft = fft(x, M);

    [peak, pos] = max(x_fft);
    %freqs = [0:Fs/(length(x_fft)-1):Fs];
    %plot(freqs, abs(x_fft));
    %omega_0_est(j) = freqs(pos);
    omega_0_est(j) = 2*pi*(pos)/(M*T);
    
    %%%%y = fftshift(x_fft);
    %%%%[my, mpos] = max(y);
    
    %%%%angles = angle(x_fft);
    %%%%angle_est(j) = angles(pos);
    angle_est(j) = angle(exp(-1i*omega_0_est(j)*best_n0(N)*T)*x_fft(pos-1));
    if mod(j, steps/100) == 0
        fprintf('%i%%\n',100* j/steps)
    end
end
clc;
err = mean(omega_0 - omega_0_est);
omega_0_var = var(omega_0_est - omega_0);
fprintf('Estimated angle: %f deg\n', rad2deg(mean(angle_est)))
fprintf('Estimated omega_0: %f\nError variance: %f | CRLB: %f\n', mean(omega_0_est), omega_0_var, best_freq(noise_std, A, T, N));

function snr = SNR(A, sigma)
    snr = A^2/(2*sigma^2);
end

function omega = best_freq(sigma, A, T, N) 
omega = 12*sigma^2 /(A^2 * T^2 * N * (N^2 - 1));
end
function phase = get_phase(sigma, n0, A, N) 
       P = make_P(N);
       Q = make_Q(N);
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