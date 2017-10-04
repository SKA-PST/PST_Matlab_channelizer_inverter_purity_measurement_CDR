
%% GLOBAL PARAMETERS

% PFB parameters
N = 256;  % number of output channels
OS_Nu = 4;  % OS factor numerator
OS_De = 3;  % OS factor denominator
PFBarms = (N*OS_De)/OS_Nu;  % number of PFB filter arms

% Length of forward FFT to process fine channels
fftLength = 2^8;

% Length of test vector block
blockLength = fftLength*N*OS_De/OS_Nu;

equaliseRipple = 1;  % 1 to equalise PFB ripple, 0 to not

% Metrics
amplitudeTF = zeros(blockLength,1);
powerTFdB = zeros(blockLength,1);
errorTFdB = zeros(blockLength,1);
maxSpuriousPower = zeros(blockLength,1);
totalSpuriousPower = zeros(blockLength,1);
averageSpuriousPower = complex(zeros(blockLength,1));


%% DESIGN PFB PROTOTYPE FILTER

%% PFB SETUP

% PFB prototype filter coefficients

% 256 channels, 4/3 OS, 4096 taps (16 per channel), 18-bit quantization
%load('CBFLowFilterCoefficients.mat');
% h = double(Num)./(32711900.0);  % scaling to give unity output amplitude

% 256 channels, 4/3 OS, 3073 taps (12 per channel + 1), 16-bit quantization
load('PSTLowFilterCoefficients.mat');

Ntaps = length(h);
    
% set up multiplicative deripple gain vector
if(equaliseRipple)
    % Produce a sampled version of the Transfer Function for later equalisation
    % - length should be N times the half-channel width (where width is FFTlength/OS_factor)
    passbandLength = ((fftLength/2)*OS_De)/OS_Nu;
    [H0,W] = freqz (h, 1, N*passbandLength);

    % use just the baseband passband section of transfer function
    % - apply to both halves of channel
    deripple = ones(passbandLength+1,1)./abs(H0(1:passbandLength+1,1));
end;


%% MAIN LOOP OVER ALL TEST FREQUENCIES

for freqBinOffset = 1 : blockLength/128
    
    % Print loop number
    fprintf('\nFrequency %d of %d\n',freqBinOffset,blockLength);
    
    % GENERATE TEST VECTOR (input to PFB)
    % amplitude 1.0
    disp('generating test vector');
    phase = rand*2.0*pi;  % randomise the phase of each sinusoid so the later summation of spectra is representative of white noise
    %phase = 0.0;          % same phase for all sinusoids -> equivalent to impulse input
    blockSig = complexSinusoidSingleBin(blockLength,freqBinOffset,phase);
    
	% generate double block length: will use middle half, after PFB primed
    inputSig = [blockSig; blockSig];  % can be concatenated cleanly because tone in bin centre

%     figure;
%     subplot(211); plot((1+blockLength/2:3*blockLength/2),real(inputSig(1+blockLength/2:3*blockLength/2,1))); box on; grid on; title('sig Real'); 
%     subplot(212); plot((1+blockLength/2:3*blockLength/2),imag(inputSig(1+blockLength/2:3*blockLength/2,1))); box on; grid on; title('sig Imag'); xlabel('time');

    usedSig = complex(zeros(blockLength,1));
    usedSig(1:blockLength,1) = inputSig(1+blockLength/2:3*blockLength/2,1);
    spectrum = fft(usedSig)./blockLength;

%     figure;
%     subplot(211); plot((1:blockLength),abs(spectrum(1:blockLength,1))); box on; grid on; title('sig Mag'); 
%     subplot(212); plot((1:blockLength),phase(spectrum(1:blockLength,1))); box on; grid on; title('sig Phase'); xlabel('frequency');


    %% PFB Channelize - one block

    disp('channelising test vector');

    % derive PFB output fine channel vector length
    NFineOut = floor(blockLength*OS_Nu/OS_De/N);

    % initialize output array
    output = complex(zeros(2*NFineOut,N));

    % channelize using CSIRO PFB analysis function
    output = polyphase_analysis(transpose(inputSig),h,N,PFBarms);
    
    % offset to start of plot/fft
    inputOffset = fftLength/2;

    disp('calculating fine channel spectra');

    % analyse middle half: length = fftLength
    samples = output(inputOffset+1:inputOffset+fftLength,:);
    spectra = fft(samples)./fftLength;
    spectra = fftshift(spectra,1);
    spectra = fftshift(spectra,2);


    %% DISCARD OVERSAMPLED PORTIONS & OPTIONALLY EQUALISE RIPPLE

    fprintf('discarding oversampled portions');
    if (equaliseRipple)
        fprintf(' and equalising ripple');
    end
    fprintf('\n');
    FN = complex(zeros(fftLength*OS_De/OS_Nu,N));
    for chan = 1:N
        % Keep only the pass-band
        discard = (1.0 - (OS_De/OS_Nu))/2.0;
        FN(:,chan) = spectra(round(discard*fftLength)+1:round((1.0-discard)*fftLength),chan);

        if(equaliseRipple)
            for ii = 1:passbandLength
                FN(ii,chan) = FN(ii,chan)*deripple(passbandLength-ii+2);
                FN(passbandLength+ii,chan) = FN(passbandLength+ii,chan)*deripple(ii);
            end;
        end;
    end;


    %% Combine chunks & back-transform
    
    fprintf('combining channels and back transforming\n');
    FN_width = fftLength*OS_De/OS_Nu;
    FFFF = FN(FN_width/2+1:FN_width,1); % upper half of chan 1 is first part of FFFF
    for chan = 2 : N
        FFFF = [FFFF; FN(:,chan)];
    end;
    FFFF = [FFFF; FN(1:FN_width/2,1)]; % lower half of chan 1 is last part of FFFF

    len = length(FFFF);
     
  	% back transform
    z1 = ifft(fftshift(FFFF))./(OS_Nu/OS_De);  % re-scale by OS factor

    if 0
        figure;
        subplot(211); plot((1:len),abs(FFFF)); box on; grid on; title('FFFF Mag'); 
        subplot(212); plot((1:len),angle(FFFF)); box on; grid on; title('FFFF Phase'); xlabel('time');
        figure;
        subplot(211); plot((1:len),real(z1(1:len))); box on; grid on; title('z1 Real'); 
        subplot(212); plot((1:len),imag(z1(1:len))); box on; grid on; title('z1 Imag'); xlabel('time');
        
        figure;
        plot((1:len),20.0*log10(abs(FFFF(1:len))+0.00000000000001)); box on; grid on; title('Output Power By Frequency'); xlabel('Frequency Bin'); ylabel('Power (dB)');
        
        fprintf('\nPress any key for next frequency...\n');
        pause;
    end;
    
    amplitudeTF(freqBinOffset,1) = abs(FFFF(freqBinOffset,1));
    FFFF(freqBinOffset,1) = 0.0;
    maxSpuriousPower(freqBinOffset,1) = (max(abs(FFFF)))^2;
    for ii = 1:blockLength
        totalSpuriousPower(freqBinOffset,1) = totalSpuriousPower(freqBinOffset,1) + (abs(FFFF(ii,1)))^2;
    end;
    averageSpuriousPower = averageSpuriousPower + abs(FFFF).^2;
    
end;

averageSpuriousPower = averageSpuriousPower./blockLength;  % averaged over number of frequencies tested

fprintf('\nMax Average Spurious Power = %3.2f dB\n',10.0*log10(max(averageSpuriousPower)));
fprintf('\nSNR = %3.2f dB\n',10.0*log10(1.0/sum(averageSpuriousPower)));  % total signal power = 1 (per frequency)

figure;
plot((1:blockLength),amplitudeTF(1:blockLength,1)); box on; grid on; axis([1 blockLength 0.999 1.001]); title('Amplitude vs Frequency'); xlabel('Frequency Bin'); ylabel('Output Amplitude');

powerTFdB = 20.0*log10(amplitudeTF);
figure;
plot((1:blockLength),powerTFdB(1:blockLength,1)); box on; grid on; axis([1 blockLength -0.1 0.1]); title('Normalised Power vs Frequency'); xlabel('Frequency Bin'); ylabel('Output Power (dB)');

powerTF = amplitudeTF.^2;
meanPower = mean(powerTF);
errorTFdB = 10.0*(log10(abs(powerTF - meanPower)) - log10(meanPower));
figure;
plot((1:blockLength),errorTFdB(1:blockLength,1)); box on; grid on; axis([1 blockLength -70 -30]); title('Power Error vs Frequency'); xlabel('Frequency Bin'); ylabel('Power (dB)');

maxSpuriousPower = 10.0*log10(maxSpuriousPower);
figure;
plot((1:blockLength),maxSpuriousPower(1:blockLength,1)); box on; grid on; axis([1 blockLength -100 -50]); title('Max Spurious Power'); xlabel('Frequency Bin'); ylabel('Max Spurious Power Component (dB)');

totalSpuriousPower = 10.0*log10(totalSpuriousPower);
figure;
plot((1:blockLength),totalSpuriousPower(1:blockLength,1)); box on; grid on; axis([1 blockLength -100 -50]); title('Total Spurious Power'); xlabel('Frequency Bin'); ylabel('Total Spurious Power (dB)');

averageSpuriousPower = 10.0*log10(averageSpuriousPower);
figure;
plot((1:blockLength),averageSpuriousPower(1:blockLength,1)); box on; grid on; axis([1 blockLength -150 -50]); title('Average Spurious Power'); xlabel('Frequency Bin'); ylabel('Average Spurious Power (dB)');


fprintf('\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
