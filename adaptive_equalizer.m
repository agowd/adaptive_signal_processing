clear close all;
%% build adaptive equalizing filter

% training sequence
N = 1000;
s = randi([0 1], 1, N) * 2 - 1;

% perform channel actions + noise
snr = 30;
h = [0.3 1 0.7 0.3 0.2];
% plotting(h);
hs = conv(h,s);
x = awgn(hs, snr);

% adaptive equalization
M = 20;
mu = 0.01;

g = zeros(1,M+1); % equalizer
e = zeros(1,length(x)); % track error

% perform LMS
num_delay = floor(M/2) + floor((length(h)-1)/2);
s_d = [zeros(1,num_delay) s]; % delayed training sequence

for i = M+1 : length(x)
    flipped_x = fliplr(x(i-M:i));
    e(i) = s_d(i) - g * flipped_x.';
    g = g + mu * e(i) * flipped_x;
end

%% training plots

output = conv(x,g);
dec = sign(output);
num_rem = floor(M/2) + floor((length(h)-1)/2);
dec_rem = dec(num_rem+1:end-num_rem);

figure;
subplot(4,1,1);
plot(s);
title('Training Sequence');
subplot(4,1,2);
plot(x);
title('x[n], signal after channel and with noise');
subplot(4,1,3);
plot(output);
title('Output from passing x[n] through Trained Equalizer');
subplot(4,1,4);
plot(dec_rem ~= s);
title('Errors between Decisions of Equalized x[n] and s[n]');
xlabel('n');

%% plot error convergence
figure;
plot(e.^2);
xlabel('Iteration');
title('Squared Error Convergence of Adaptive Equalization Filter');

%% plots of channel
plotting(h);
figure;
zplane(h);
title('Pole-Zero Plot of Channel')

plotting(g);

plotting(conv(h,g));

%% testing filter – generate test sequence
test_length = 10000;
s_test = randi([0 1], 1, test_length) * 2 - 1;
hs_test = conv(h,s_test);
x_test = awgn(hs_test, snr);

% testing filter – filter test sequence
filtered_x = conv(x_test, g);
num_delay = floor(M/2) + floor(length(h)/2);
remove_extras = filtered_x(num_delay+1:end-num_delay);

x_test_rem = x_test(floor(length(h)/2)+1:end-floor(length(h)/2));
equalized_output_decisions = sign(remove_extras);
unequalized_output_decisions = sign(x_test_rem);

% compare filtered results to actual
figure;
subplot(1,2,1);
plot((x_test_rem - s_test).^2);
xlabel('n');
title('Squared Error Between Un-Equalized and Signal');
ylim([0 25]);

subplot(1,2,2);
plot((remove_extras - s_test).^2);
xlabel('n');
title('Squared Error Between Equalized and Signal');
ylim([0 25]);

total_incorrect_unequalized_output_decision = sum(unequalized_output_decisions ~= s_test);
total_incorrect_equalized_output_decision = sum(equalized_output_decisions ~= s_test);

rate_uneq = total_incorrect_unequalized_output_decision/test_length * 100
rate_eq = total_incorrect_equalized_output_decision/test_length * 100

times = (1-rate_eq/100) / (1 - rate_uneq/100)

%% function for plotting + freq
function plotting(x)
    n = 512;
    [freq,w] = freqz(x,n);

    figure;
    subplot(2,4,[1,2,5,6]);
    stem(x);
    xlabel('n');
    title('Impulse Response');
    subplot(2,4,[3,4]);
    plot(w/pi, abs(freq));
    xlabel('Angular Frequency (\times\pi rad/sample)');
    title('Magnitude Response')
    subplot(2,4,[7,8]);
    plot(w/pi, angle(freq));
    xlabel('Angular Frequency (\times\pi rad/sample)');
    title('Phase Response')
end