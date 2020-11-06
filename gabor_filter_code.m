clc; clear all; close all;

% Load and plot the The Messiah segment

load handel
v = y'/2;
plot((1:length(v))/Fs,v);
xlabel("Time [sec]");
ylabel("Amplitude");
title("Handel's Messiah");
axis([0 9 0 1])
%% Plot with a 'window'
t_mesh = (1:length(v))/Fs;
hold on
plot(t_mesh, exp(-(t_mesh - 4).^2), 'LineWidth', 2)
legend('Spatial Signal', "G\'abor Window", 'Interpreter', 'latex')

%% Playback

p8 = audioplayer(v,Fs);
playblocking(p8);

%% Build the time and frequency domains

L = length(v)/Fs; n = length(v);
t2 = linspace(0, L, n+1); t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k);
v = v(1:end - 1);
v_t = fft(v);

%% Plot in freq. domain of original file
close all
plot(ks,abs(fftshift(v_t))/max(abs(v_t)),'k'); %axis([-50 50 0 1])
%set(gca,'Fontsize',[14])
xlabel('frequency (\omega)'), ylabel('Amplitude')

%% Gabor Window Iteration/Slide
close all
%figure(4)
tslide=0:0.1:L;

messiah_spec_1 = zeros(length(tslide), length(t)-1);
messiah_spec_2 = zeros(length(tslide), length(t)-1);

% Two different bands to compare
w1 = .1;
w2 = 100;

for j=1:length(tslide)
    g=exp(-w1*(t(1:end-1)-tslide(j)).^2); % Gabor 
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_spec_1(j, :) = abs(fftshift(v_g_t)); 
end
for j=1:length(tslide)
    g=exp(-w2*(t(1:end-1)-tslide(j)).^2); % Gabor 
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_spec_2(j,:) = abs(fftshift(v_g_t)); 
end

figure(5)
subplot(2,1,1)
pcolor(tslide,ks,messiah_spec_1.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('frequency (\omega)')
title('Wide Gabor Filter');

subplot(2,1,2)
pcolor(tslide,ks,messiah_spec_2.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('frequency (\omega)')
title('Narrow Gabor Filter')

%% Under versus oversampling
close all
%figure(4)
tslide1=0:0.1:L;
tslide2=0:1:L;
messiah_spec_1 = zeros(length(tslide1), length(t)-1);
messiah_spec_2 = zeros(length(tslide2), length(t)-1);


w = 100;

for j=1:length(tslide1)
    g=exp(-w*(t(1:end-1)-tslide1(j)).^2); % Gabor 
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_spec_1(j, :) = abs(fftshift(v_g_t)); 
end
for j=1:length(tslide2)
    g=exp(-w*(t(1:end-1)-tslide2(j)).^2); % Gabor 
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_spec_2(j,:) = abs(fftshift(v_g_t)); 
end

subplot(2,1,1)
pcolor(tslide1,ks,messiah_spec_1.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('frequency (\omega)')
title('Small (t=.1) window steps');

subplot(2,1,2)
pcolor(tslide2,ks,messiah_spec_2.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('frequency (\omega)')
title('Large (t=2) window steps')

%% Plot the different filters
close all

x_mesh = linspace(-3, 3, 100);
y_g = exp(-20*(x_mesh(1:end)).^2);
y_m = (1 - x_mesh.^2).*exp(-10*(x_mesh).^2);
t_shift = t(1:end -1)- 2;
y_sh = t_shift < 1/4 & t_shift > -1/4;

subplot(3,1,1)
plot(x_mesh, y_g)
axis([-3 3 -.2 1.2])
title("Gaussian Filter")

subplot(3,1,2)
plot(x_mesh, y_m)
axis([-3 3 -1 1.2])
title("Mexican Hat Filter")


subplot(3,1,3)
plot(t(1:end-1), y_sh)
axis([0 4 -.2 1.2])
title("Shannon Filter")


%% Plot Spectorgrams with Multiple types of filters
close all

tslide=0:0.05:L;
w = 100;

% Spatial resolution with the shannon filter is pretty outstanding
% with this filter width - the doubled 'hallelujah's are pretty clear
% in the spectrogram.
w_sh = .1;
messiah_sp_gauss = zeros(length(tslide), length(t)-1);;
messiah_sp_mex = zeros(length(tslide), length(t)-1);;
messiah_sp_shan = zeros(length(tslide), length(t)-1);;

for j=1:length(tslide)
    g=exp(-w*(t(1:end-1)-tslide(j)).^2);
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_sp_gauss(j,:)=abs(fftshift(v_g_t)); 
end

for j=1:length(tslide)
    g=(1 - (t(1:end -1) - tslide(j)).^2).*exp(-w*(t(1:end-1)-tslide(j)).^2); 
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_sp_mex(j,:)=abs(fftshift(v_g_t)); 
end

for j=1:length(tslide)
    
    % Shift the t-mesh by the slide value, then use this
    % to return a vector of booleans depending on whether or not 
    % the value in the t-mesh is greater than the -width/2 AND less 
    % than the width/2. This centers a square filter at tslide(j) or tau.
    t_shift = t(1:end -1) - tslide(j);
    g = t_shift < w_sh/2 & t_shift > -w_sh/2;
    v_g= g.*v; 
    v_g_t=fft(v_g); 
    messiah_sp_shan(j,:)=abs(fftshift(v_g_t)); 
end

%% Plot the data using pcolor
close all

subplot(3,1,1)
pcolor(tslide,ks,messiah_sp_gauss.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('\omega')
title('Gaussian Filter');

subplot(3, 1, 2)
pcolor(tslide,ks,messiah_sp_mex.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('\omega')
title('Mexican Hat Filter')

subplot(3,1,3)
pcolor(tslide,ks,messiah_sp_shan.'), 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
xlabel('time(seconds)'), ylabel('\omega')
title('Shannon Filter')

%% Part II

clc; clear all; close all;

% tr_piano=16; % record time in seconds
% y_p=audioread('music1.wav'); Fs=length(y_p)/tr_piano;
% plot((1:length(y_p))/Fs,y_p);
% xlabel('Time [sec]'); ylabel('Amplitude');
% title('Mary had a little lamb (piano)'); drawnow
% p8 = audioplayer(y_p,Fs); playblocking(p8);
% figure(2)
tr_rec=14; % record time in seconds
y_p= gpuArray(audioread('music2.wav')); Fs=length(y_p)/tr_rec;
plot((1:length(y_p))/Fs,y_p);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
% p8 = audioplayer(y_p,Fs); playblocking(p8);

%%
L = 14; n = length(y_p);
t2 = linspace(0, L, n+1); t = t2(1:n);
 % k = (2*pi/L)*[0:n/2-1 -n/2:-1]; ks = fftshift(k);
k = (1/L)*[0:(n/2-1) -n/2:-1];
ks = fftshift(k);

%% Plot the Piano in Frequency Domain
y_p_t = fft(y_p);
plot(ks, -abs(fftshift(y_p_t))/max(y_p_t), 'LineWidth', 1)
axis([750 2200 -.1 1.1])
xticks([700 816 910 1034 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200])
xticklabels({'700', '816', '910', '1034', '1100', '1200', '1300', '1400',  '1500', '1600', '1700', '1800', '1900', '2000', '2100', '2200'})
xlabel('Frequency (Hz)'), ylabel('Amplitude')
title('Recorder - Frequency Domain');
%% Check width of gabor filter
w = .1;
plot(t,exp(-w*(t(1:end)-8).^2))
 
%% Gabor Window Iteration/Slide
close all; clc;
%figure(4)
tslide=0:0.05:L;
piano_spec = zeros(length(tslide), length(t));
w = 100;
y_p = y_p';


for j=1:length(tslide)
    g=exp(-w*(t(1:end)-tslide(j)).^2); % Gabor 
    y_g= g.*y_p; 
    y_g_t=fft(y_g); 
    piano_spec(j,:)= abs(fftshift(y_g_t)); 
end

% Range of frequency to plot, based on observation of the frequency domain
% plot.
idx = find(ks > 750 & ks < 1150);

pcolor(tslide,ks(idx),piano_spec(:,idx).') 
shading interp 
set(gca,'Fontsize',[14]) 
colormap(hot)
%axis([0 16 -.2 .2])
xlabel('time(seconds)'), ylabel('frequency (Hz)')
title('Narrow Gabor Filter');

%%
y_p_t_s = fftshift(y_p_t);
% [B, I] = maxk(y_p_t_s, 3);
% This wound up being annoying because the top three max freq were
% all around the first (loudest) note.

% Gaussian filters to remove all overtones. The k1/2/3 are
% simply found by examining a frequency domain plot of the sample.

% Piano Peaks
% k1 = 255;
% k2 = 287;
% k3 = 318;

% Recorder Peaks
k1 = 816;
k2 = 910;
k3 = 1034;

w = .1;

f1 = exp( -w *(ks - k1).^2);
f2 = exp( -w *(ks - k2).^2);
f3 = exp( -w *(ks - k3).^2);

y_p_filtered = (y_p_t_s.* (f1 + f2 + f3)');
%% Plot the filters To confirm width
% Note the filter width matters, since a very narrow filter
% will drastically reduce the amplitude and produce less contrast
% in the pcolor plot with sound versus silence.
plot(ks, f1 + f2 + f3)
axis([750 1200 -.1 1.1])
xlabel('time(seconds)'), ylabel('frequency (Hz)')
title('Low-Pass Filter');
%% Plot the filtered frequency domain

close all
plot(ks, -abs(y_p_filtered)/max(y_p_filtered))
axis([750 1100 -.1 1.1])
xlabel('Frequency(Hz)'), ylabel('Amplitude')
title('Filtered Recorder - Freq. Domain');
%% Now a spectrogram of the filtered data
close all; clc;
%figure(4)
tslide=0:0.05:L;
filtered_spec = zeros(length(tslide), length(t));


% Two different bands to compare
w = .1;
w_sh = .1;

y_p_f = ifft(ifftshift(y_p_filtered));


for j=1:length(tslide)
    g=exp(-w*(t(1:end)-tslide(j)).^2); % Gabor 
    
    % Shannon Filter (gaussian seems to work better)
    % t_shift = t(1:end) - tslide(j);
    % g = t_shift < w_sh/2 & t_shift > -w_sh/2;
    
    y_g= g'.*y_p_f; 
    y_g_t=fft(y_g); 
    filtered_spec(j,:) = abs(fftshift(y_g_t));
end

%% Plot using pcolor 

pcolor(tslide,ks(idx),filtered_spec(:,idx).') 
shading interp 
set(gca,'Fontsize',[14], 'YLim', [750 1150], 'YTick', [750:20:1150]) 
% set(gca,'Fontsize',[14]) 
colormap(hot)
%axis([0 16 -.2 .2])
xlabel('time(seconds)'), ylabel('frequency (Hz)')
title('Filtered Data, Narrow Gabor Filter');

% Results are perfect. The Piano is out of tune, ranging from around 5
% to 10 Hz flat.

%% Play the filtered data

% Note the playback sounds like pure tones, just as we would expect
p8 = audioplayer(y_p_f,Fs); playblocking(p8);
