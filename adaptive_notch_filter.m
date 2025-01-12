clear, close all;

%% part 1
% generate desired signal + sinusoidal interference
t = 1;
fs = 1000;
n = linspace(0, fs * t, fs*t+1);

% sum of sinusoids
desired = sin(2 * pi * (93.75/fs) .* n);

int_freq = 187.5;
interference = 1.15*sin(2 * pi * (int_freq/fs) .* n);

x = desired+interference;
y = zeros(size(x));

% perform adaptive filtering
a = 0;
mu = 0.001;
r = 0.98;

track_a = zeros(size(x));
error = zeros(size(x));

for i = n(3:end)+1
    e = x(i) + a * x(i-1) + x(i-2);
    y(i) = e - r * a * y(i-1) - r^2 * y(i-2);
    a = a - mu * y(i) * x(i-1);
    error(i) = y(i)-desired(i);
    if a < -2 || a > 2
        a = 0;
    end
    track_a(i) = a;
end

plotting(desired);
sgtitle('Desired Signal s[n]');
plotting(interference);
sgtitle('Interfering Sinusoid w[n]');
plotting(x);
sgtitle('Signal with Interference x[n]');
plotting(y);
sgtitle('Output y[n]')
plotting(y(700:end));
sgtitle('Output y[n] (after convergence)')


figure;
plot(track_a);
xlabel('n');
ylabel('a');
title('a as a function of n')

figure;
plot(error.^2);
xlabel('n');
ylabel('(y[n]-s[n])^2')
title('Squared Error for Each Iteration')

b = [1 track_a(end) 1];
a = [1 r * track_a(end) r^2];

[H, w] = freqz(b, a, 'whole');
figure;
plot(w-pi, fftshift(abs(H)));
xlabel('\omega')
ylabel('Magnitude')
title('Notch Filter Magnitude Response')

op = filter(b, a, x);
plotting(op);
sgtitle('Filtered Using Coefficients and filter() in MATLAB')
%% part 2

t = 1;
fs = 1000;
n = linspace(0, fs * t, fs*t+1);

% sum of sinusoids
desired = sin(2 * pi * (93.75/fs) .* n);

% interference = chirp(t, 0, t(end), 187.5);
interference = 8*sin(pi .* 100 .* (n./fs).^2 + 2 .* pi .* (187.5/fs) .* n);
% interference = 2*sin(2 * pi * (187.5/fs) .* n);

x = desired+interference;

y = zeros(size(x));

% perform adaptive filtering
a = 0;

mu = 0.01;
r = 0.85;

track_a = ones(size(x)) * a;
error = zeros(size(x));


for i = n(3:end)+1
    e = x(i) + a * x(i-1) + x(i-2);
    y(i) = e - r * a * y(i-1) - r^2 * y(i-2);
    a = a - mu * y(i) * x(i-1);
    error(i) = y(i)-desired(i);
    if a < -2 || a > 2
        a = 0;
    end
    track_a(i) = a;
end




plotting(desired);
sgtitle('Desired Signal s[n]');
plotting(interference);
sgtitle('Interfering Sinusoid w[n]');
plotting(x);
sgtitle('Signal with Interference x[n]');
plotting(y);
sgtitle('Output y[n]')

figure;
plot(track_a);
xlabel('n');
ylabel('a');
title('a as a function of n')

figure;
plot(error.^2);
xlabel('n');
ylabel('(y[n]-s[n])^2')
title('Squared Error for Each Iteration')

%% part 3
% generate desired signal + sinusoidal interference
t = 1;
fs = 1000;
n = linspace(0, fs * t, fs*t+1);

% sum of sinusoids
desired = sin(2 * pi * (125/fs) .* n);

int_freq_1 = 1000 / 6;
interference_1 = 3*sin(2 * pi * (int_freq_1/fs) .* n);
int_freq_2 = 1000 / 12;
interference_2 = 3*sin(2 * pi * (int_freq_2/fs) .* n);


x = desired + interference_1 + interference_2;
y_1 = zeros(size(x));
y_2 = zeros(size(x));

% perform adaptive filtering
a_1 = 0;
a_2 = 0;

mu = 0.0008;
r = 0.85;

track_a1 = zeros(size(x));
track_a2 = zeros(size(x));

% first get rid of I1
for i = n(3:end)+1
    e = x(i) + a_1 * x(i-1) + x(i-2);
    y_1(i) = e - r * a_1 * y_1(i-1) - r^2 * y_1(i-2);
    a_1 = a_1 - mu * y_1(i) * x(i-1);
    if a_1 < -2 || a_1 > 2
        a_1 = 0;
    end
    track_a1(i) = a_1;
end

% filter out I1 to be left w D + I2
b1 = [1 track_a1(end) 1];
a1 = [1 r * track_a1(end) r^2];


op1 = filter(b1, a1, x);

% find a2 from D + I2
for i = n(3:end)+1
    e = op1(i) + a_2 * op1(i-1) + op1(i-2);
    y_2(i) = e - r * a_2 * y_2(i-1) - r^2 * y_2(i-2);
    a_2 = a_2 - mu * y_2(i) * op1(i-1);
    if a_2 < -2 || a_2 > 2
        a_2 = 0;
    end
    track_a2(i) = a_2;
end

% filter out I2 from D + I1 + I2 to be left with D + I1
b2 = [1 track_a2(end) 1];
a2 = [1 r * track_a2(end) r^2];

[H2, w] = freqz(b2, a2, 'whole');
% figure;
% plot(w-pi, fftshift(abs(H2)));
% xlabel('\omega')
% ylabel('Magnitude')
% title('Notch Filter Magnitude Response')

op2 = filter(b2, a2, x);

track_a11 = zeros(size(x));

% find a1 from D + I1
for i = n(3:end)+1
    e = op2(i) + a_1 * op2(i-1) + op2(i-2);
    y_1(i) = e - r * a_1 * y_1(i-1) - r^2 * y_1(i-2);
    a_1 = a_1 - mu * y_1(i) * op2(i-1);
    if a_1 < -2 || a_1 > 2
        a_1 = 0;
    end
    track_a11(i) = a_1;
end


b11 = [1 track_a11(end) 1];
a11 = [1 r * track_a11(end) r^2];

[H11, w] = freqz(b11, a11, 'whole');

op_ = filter(b11, a11, x);
op_ = filter(b2, a2, op_);

% bfr = conv(b11,b2);
% afr = conv(a11, a2);
% op_ = filter(bfr, afr, x);

plotting(desired);
sgtitle('Desired Signal s[n]');
plotting(interference_1);
sgtitle('Interfering Sinusoid 1 w_1[n]');
plotting(interference_2);
sgtitle('Interfering Sinusoid 2 w_2[n]');
plotting(x);
sgtitle('Signal with Interference x[n]');
plotting(op_);
sgtitle('Filtered Output using filter() in MATLAB');

figure;
plot(track_a11);
xlabel('n');
ylabel('a');
title('a_1 as a function of n')

figure;
plot(track_a2);
xlabel('n');
ylabel('a');
title('a_2 as a function of n')

figure;
plot((w-pi),fftshift(abs(H11).*abs(H2)));
xlabel('\omega')
ylabel('Magnitude')
title('Notch Filter Magnitude Response')

figure;
plot((op_-desired).^2);
xlabel('n');
ylabel('(y[n]-s[n])^2')
title('Squared Error for Each Sample')

%% function for plotting + freq
function [Y, w, gd, w1] = plotting(x)
    L = 1001;
    Y = fft(x,L);
    w = linspace(-1*pi,pi,L);
    n = 1:length(x);
    figure;
    subplot(3,1,1);
    plot(n,x);
    xlabel('n');
    title('Time Domain')
    subplot(3,1,2);
    plot(w/pi,fftshift(abs(Y)));
    xlabel('\omega (pi * rad/sample)');
    title('Magnitude');
    subplot(3,1,3);
    plot(w/pi,fftshift(angle(Y)));
    xlabel('\omega (pi * rad/sample)');
    title('Phase');
end