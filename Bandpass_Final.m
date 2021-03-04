close all;
clear;
clc

%prescribed parameters A,B, and C are calculated by the user input

indexIn = inputdlg ( ' Enter Index Number' , 'Index Number' , 1 )  ;
indexNo = str2double (indexIn{1}(1:end-1)); 
C = mod(indexNo , 10 ) ;
B = mod( floor(indexNo/10) ,10 ) ;
A = mod( floor(indexNo/100) ,10 ) ;

%parameters for the Bandpass filter

A_tilda_p=0.03+0.01*A;%Maximum passband ripple-:0.0900
A_tilda_a=45+B;%Minimum stopband attenuation-:45
wp1=C*100+300;%Lower passband edge-:1200
wp2=C*100+700;%Upper passband edge-:1600
wa1=C*100+150;%Lower stopband edge-:1050
wa2=C*100+800;%Upper stopband edge-:1700
ws=2*(C*100+1200);%Sampling frequency-:4200

% Derived Specifications

bt1 = wp1-wa1;  %lower transition width
bt2 = wa2-wp2;  % upper transisiton width
bt = min(bt1,bt2); %critical transition width
wc1 = wp1-bt/2; % lower cutoff frequency
wc2 = wp2+bt/2; % upper cutoff frequency
T = 2*pi/ws; % sampling period

% Kaiser Window Parameters

deltaP = (10^(0.05*A_tilda_p) - 1)/ (10^(0.05*A_tilda_p) + 1); % calculating delta
deltaA = 10^(-0.05*A_tilda_a);
delta = min(deltaP,deltaA);

Aa = -20*log10(delta);  % Actual stopband attenuation

if Aa<=21               % Calculating alpha 
    alpha = 0;
elseif Aa>21 && Aa<= 50
    alpha = 0.5842*(Aa-21)^0.4 + 0.07886*(Aa-21);
else
    alpha = 0.1102*(Aa-8.7);
end

if Aa <= 21             % Calculating D
    D = 0.9222;
else
    D = (Aa-7.95)/14.36;
end

N = ceil(ws*D/bt +1); % order of the filter
if mod(N,2) == 0
    N = N+1;
end

n = -(N-1)/2:1:(N-1)/2;   % length of the filter

beta = alpha*sqrt(1-(2*n/(N-1)).^2);

% Generating I(alpha)

bessellimit = 50;

Ialpha = 1;
for k = 1:bessellimit
    termk = (1/factorial(k)*(alpha/2).^k).^2;
    Ialpha = Ialpha + termk;
end

% Generating I(beta) %%

Ibeta = 1;
for k = 1:bessellimit
    termk = (1/factorial(k)*(beta/2).^k).^2;
    Ibeta = Ibeta + termk;
end

% Obtaining Kaiser Window 

wknt = Ibeta/Ialpha;

figure
stem(n,wknt)
xlabel('n')
ylabel('Amplitude')
title('Kaiser Window - Time Domain');

% Generating Impulse Response 

nleft = -(N-1)/2:-1;
hntleft = 1./(nleft*pi).*(sin(wc2*nleft*T)-sin(wc1*nleft*T));

nright = 1:(N-1)/2;
hntright = 1./(nright*pi).*(sin(wc2*nright*T)-sin(wc1*nright*T));

hnt0 = 2/ws*(wc2-wc1);

hnt = [hntleft,hnt0,hntright];
figure
stem(n,hnt)
xlabel('n')
ylabel('Amplitude')
title(strcat('Filter Response - Rectangular window - Time Domain'));

figure
[hi_f,wi_f] = freqz(hnt);
wi_f = wi_f/T;
hi_f = 20*log10(abs(hi_f));
plot(wi_f,hi_f)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat('Filter Response - Rectangular Window - Frequency Domain'));

% Applying the window to the filter


filter = hnt.*wknt;
figure
stem(n,filter)
xlabel('n')
ylabel('Amplitude')
title(strcat('Filter Response - Kaiser Window - Time Domain'));

figure
[h,w] = freqz(filter);
w = w/T;
h = 20*log10(abs(h));
plot(w,h)
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title(strcat('Filter Response - Kaiser Window - Frequency Domain'));

% Plotting the Passband %%

figure
start = round(length(w)/(ws/2)*wc1);
finish = round((length(w)/(ws/2)*wc2));
wpass = w(start:finish);
hpass = (h(start:finish));
plot(wpass,hpass)
axis([-inf, inf, -0.1, 0.1]);
xlabel('Frequency (rad/s)')
ylabel('Magnitude (dB)')
title('Passband - Frequency Domain');



% Input signal generation %%

w1 = wc1/2;                    %% component frequencies of the input
w2 = wc1 + (wc2-wc1)/2;
w3 = wc2 + (ws/2-wc2)/2;
wi = [w1,w2,w3];

n1 = 0:1:500;                  %% x axis for discrete signal
n2 = 0:0.1:500;                %% x axis for envelope

xnt = 0;     
xdash =0;

for each = wi                        
    xnt = xnt+ sin(each.*n1.*T);     %% generate discrete signal
    xdash = xdash+sin(each.*n2.*T);  %% generate envelope
end


% Using DFT to check the filtering %%

% Filtering using frequency domain multiplication - see function getfiltered.m %%

Npoint = length(xnt) + length(filter) - 1;  % length for fft in x dimension
xfft = fft(xnt,Npoint);
filterfft = fft(filter,Npoint);
outfft = filterfft .* xfft;
out = ifft(outfft,Npoint);

out1 = out(floor(N/2)+1:length(out)-floor(N/2)); % account for shifting delay


% Frequency domain representation of input signal before filtering

figure 
subplot(2,1,1)

Npoint = length(xnt);  
xfft = fft(xnt,Npoint);
x1 = n1/length(n1)*ws-ws/2;
plot(x1,abs(xfft))
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat(['Input signal',' ','- Frequency Domain']));

% Frequency domain representation of input signal after filtering

subplot(2,1,2)
Npoint = length(out1);  % length for fft in x dimension
xffout = fft(out1,Npoint);
x1 = n1/length(n1)*ws-ws/2;
plot(x1,abs(xffout))
xlabel('Frequency rad/s')
ylabel('Magnitude')
title(strcat('Output signal',' ','- Frequency Domain'));



% time domain representation of input signal before filtering

figure
subplot(2,1,1)
stem(n1,xnt)
xlabel('n')
ylabel('Amplitude')
title(strcat(['Input signal',' ','- Time Domain']));
hold on 
plot(n2,xdash)
legend('Input signal','Input signal envelope');



% Time domain representation of input signal after filtering

subplot(2,1,2)
stem(n1,out1)
xlabel('n')
ylabel('Amplitude')
title(strcat('Output signal',' ','- Time Domain'));

hold on
plot(n2,sin(w2.*n2.*T))
legend('Output signal','Envelope of the middle frequency component of the input signal');


