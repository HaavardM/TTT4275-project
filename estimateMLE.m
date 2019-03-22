A = 1;
phase = pi/8;
f0 = 10^5;
omega_0 = 2*pi*f0;
Fs = 10^6;
T = 1/Fs;

N = 513;
noise_var = 1;

signal = A*exp(1i*(omega_0*[1:N]*T + phase));
noise = normrnd(0, noise_var, 1, N) + 1i*normrnd(0, noise_var, 1, N);
x = signal + noise;

x_fft = fft(x);


[peak, pos] = max(abs(x_fft(2:end)));
freqs = [0:Fs/(N-1):Fs];
plot(freqs, abs(x_fft));
omega_0_est = freqs(pos);
fprintf('Estimated omega_0: %f', omega_0_est);

function snr = SNR(A, sigma)
    snr = A^2/(2*sigma^2);
end

function omega = var_freq(sigma, A, T, N) 
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
    n0 = make_no(make_P(N), N);
end
function n0 = make_n0(P, N)
    n0 = -P/N;
end