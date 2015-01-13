%Unprincipled FFT. Just powerful enough to get you in trouble.
%Wild type FFT
a = rand(1,629);
b = 0:0.1:20*pi;
c = sin(5*b);
d = a+c;
figure
subplot(1,2,1)
plot(b,d)

subplot(1,2,2)
e = fft(d)
plot(b,e)

%%
%Principled FFT. Developing it from first principles
%Create the time series successively 
figure
L = 1000; %Length of the signal. 1000 paces
X = zeros(1,L); %Start with zeros
sampling_interval = 0.1; %We sample 10 times per second. 10 Hz.
t = (1:L) * sampling_interval; %Make a time series, sampled at the frequency
ns = 10; %number of sines to add
for nn = 1:ns %Add the suitable number of sine waves to the signal
    X = X + nn * sin (nn*pi*t); %Add it all together. Each pass we add a sine wave
end
subplot(1,2,1)
plot(t,X); %Let's have a look what the time series looks like
Y = fft(X)/L; %Take the fft of the signal, normalize by length
NyLimit = (1/sampling_interval)/2; %Calculate Nyquist
F = linspace(0,1,L/2)*NyLimit;%Understanding the properties of the FFT, we don't want to plot from -Nyquist to +Nyquist, we only want to plot half of it. The other half is symmetric
subplot(1,2,2)
plot(F,abs(Y(1:L/2))); %Plotting only half of the space. Other half is mirror symmetric and completely redundant

%Power considerations. 
%Power is basically the frequency representation of the signal multiplied
%by its complex conjugate
power_of_Y = Y.*conj(Y); %Using the definition of power, we calculate it
figure
subplot(1,2,1)
plot(F,power_of_Y(1:L/2)) %Plot the power to feel the power

phi = atan(imag(Y)./real(Y)); %Determine the phase angles of the signal
subplot(1,2,2)
plot(F,phi(1:L/2)); %Plotting the series of the phase angles

%%
%Doing it for real. And imaginary.
[y_1, fs_1, nbits_1] = wavread('voice-ah.wav'); %Read in the ah sound, recognizing it is a wav file
[y_2, fs_2, nbits_2] = wavread('voice-oh.wav'); %Read in the oh sound, recognizing it is a wav file

%What does it sound like?
%Commented out due to annoying
%sound(y_1,fs_1); %Play the "ah"
%pause
%sound(y_2,fs_2); %Play the "oh"
 
%What does it look like?
figure
subplot(1,2,1)
plot(y_1)
subplot(1,2,2)
plot(y_2)

%That does not tell us anything. So let's consider the power at different
%frequencies
New_nyquist = fs_1/2; %That is our new nyquist
%Determine the length of the signal
New_F = linspace(0,1,length(y_1)/2)*New_nyquist; %Get a new basis
New_frequency_signal = fft(y_1)/length(y_1);%Don't forget to take the FFT
New_power = New_frequency_signal.*conj(New_frequency_signal); %For the first one 
figure
subplot(1,2,1)
plot(New_F,New_power(1:length(New_power)/2))
subplot(1,2,2)
%Do the same thing for the "oh"
New_F = linspace(0,1,length(y_2)/2)*New_nyquist; %Get a new basis
New_frequency_signal = fft(y_2)/length(y_2);%Don't forget to take the FFT
New_power = New_frequency_signal.*conj(New_frequency_signal); %For the first one 
plot(New_F,New_power(1:length(New_power)/2))

%%
%Spectrogram
%Why do we need it? 
%Let's load in the birdsong to appreciate why we need the spectrogram.
%A spectrogram is not needed if the signal is stationary. Alas, in biology
%almost all signals are not stationary. 
[birdsong,fs,nbits] = wavread('song1.wav');
%What does it look and sound like
figure
plot(birdsong)
sound(birdsong,fs)
figure
spectrogram(birdsong, 256, 'yaxis')

%%
%Let's first do the spectrogram of something that is also not stationary,
%but that we understand
Fs = 1000; %Sampling frequency is a thousand
t = 0:0.001:2; %0 to 2 seconds at 1 kHz
y = chirp(t,100,1,200,'quadratic'); %Start at 100 Hz, at 1 second cross 200 Hz
spectrogram(y,hamming(512),511,128,Fs,'yaxis')

