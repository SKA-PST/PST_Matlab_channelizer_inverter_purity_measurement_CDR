% Script to test oversampled PFB inversion via FFT for SKA_Low channelization
% and measure the resulting temporal leakage

% Ian Morrison
% October 2017

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
amplitudeTF = zeros(PFBarms,1);
powerTFdB = zeros(PFBarms,1);
errorTFdB = zeros(PFBarms,1);
maxSpuriousPower = zeros(PFBarms,1);
totalSpuriousPower = zeros(PFBarms,1);
averageSpuriousPower = zeros(blockLength,1);


%% PFB SETUP

% PFB prototype filter coefficients

% 256 channels, 4/3 OS, 2561 taps (10 per channel + 1), 16-bit quantization
load('PST_2561_LowFilterCoefficients.mat');

% 256 channels, 4/3 OS, 3073 taps (12 per channel + 1), 16-bit quantization
%load('PST_3073_LowFilterCoefficients.mat');

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


%% MAIN LOOP OVER ALL PFB ARMS

for impulseOffset = 1 : PFBarms
    
    % Print loop number
    fprintf('\nPFB Arm %d of %d\n',impulseOffset,PFBarms);

    % GENERATE TEST VECTOR (input to PFB)
    % impulse of amplitude 1.0
    % generate double block length: will use middle half, after PFB primed
    disp('generating test vector');
    inputSig = complex(zeros(2*blockLength,1,'single'));
    impulse_offset = blockLength + impulseOffset;  % location of impulse (first is at blockLength + 1)
    impulse_width = 1;  % number of samples width of impusle
    inputSig(impulse_offset:impulse_offset+impulse_width-1,1) = 1.0;

    % figure;
    % subplot(211); plot((1:2*blockLength),real(inputSig(1:2*blockLength,1))); box on; grid on; title('sig Real'); 
    % subplot(212); plot((1:2*blockLength),imag(inputSig(1:2*blockLength,1))); box on; grid on; title('sig Imag'); xlabel('time');

    spectrum = fft(inputSig)./(2*blockLength);

    % figure;
    % subplot(211); plot((1:2*blockLength),abs(spectrum(1:2*blockLength,1))); box on; grid on; title('sig Mag'); 
    % subplot(212); plot((1:2*blockLength),phase(spectrum(1:2*blockLength,1))); box on; grid on; title('sig Phase'); xlabel('frequency');


    %% PFB Channelize - one block

    disp('channelising test vector');

    % derive PFB output fine channel vector length
    NFineOut = floor(blockLength*OS_Nu/OS_De/N);

    % initialize output array
    output = complex(zeros(2*NFineOut,N));

    % channelize using CSIRO PFB analysis function
    output = polyphase_analysis(transpose(inputSig),h,N,PFBarms);

    % scale up by number of channels
    output = output.*N;
    
    % offset to start of plot/fft
    inputOffset = fftLength/2;

    disp('calculating fine channel spectra');

    % analyse middle half: length = fftLength
    samples = output(inputOffset+1:inputOffset+fftLength,:);
    spectra = fft(samples);
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
    
    % shift centre to common location for all impulseOffset cases
    z1 = circshift(z1, -impulseOffset);

    if 0
        figure;
        subplot(211); plot((1:len),abs(FFFF)); box on; grid on; title('FFFF Mag'); 
        subplot(212); plot((1:len),angle(FFFF)); box on; grid on; title('FFFF Phase'); xlabel('time');
        figure;
        subplot(211); plot((1:len),real(z1(1:len))); box on; grid on; title('z1 Real'); 
        subplot(212); plot((1:len),imag(z1(1:len))); box on; grid on; title('z1 Imag'); xlabel('time');

        figure;
        plot((1:len),20.0*log10(abs(z1(1:len))+0.00000000000001)); box on; grid on; title('Output Power By Sample'); xlabel('sample'); ylabel('Power (dB)');
        
        fprintf('\nPress any key for next PFB arm...\n');
        pause
    end;

    output_offset =  (OS_De*N/OS_Nu)*(fftLength-inputOffset) - floor(Ntaps/2);

    amplitudeTF(impulseOffset,1) = abs(z1(output_offset,1));
    
    z1(output_offset,1) = 0.0;
    maxSpuriousPower(impulseOffset,1) = (max(abs(z1)))^2;
      
    for ii = 1:blockLength
      totalSpuriousPower(impulseOffset,1) = totalSpuriousPower(impulseOffset,1) + (abs(z1(ii,1)))^2;
    end;
    averageSpuriousPower = averageSpuriousPower + abs(z1).^2;
        

end;

averageSpuriousPower = averageSpuriousPower./PFBarms;  % averaged over number of impulses tested

figure;
plot((1:PFBarms),amplitudeTF(1:PFBarms,1)); box on; grid on; axis([1 PFBarms 0.999 1.001]); title('Amplitude vs Impulse Location'); xlabel('Impulse Location'); ylabel('Output Amplitude');

powerTFdB = 20.0*log10(amplitudeTF);
figure;
plot((1:PFBarms),powerTFdB(1:PFBarms,1)); box on; grid on; axis([1 PFBarms -0.1 0.1]); title('Normalised Power vs Impulse Location'); xlabel('Impulse Location'); ylabel('Output Power (dB)');

powerTF = amplitudeTF.^2;
meanPower = mean(powerTF);
errorTFdB = 10.0*(log10(abs(powerTF - meanPower)) - log10(meanPower));
figure;
plot((1:PFBarms),errorTFdB(1:PFBarms,1)); box on; grid on; axis([1 PFBarms -70 -30]); title('Power Error vs Impulse Location'); xlabel('Impulse Location'); ylabel('Power (dB)');

maxSpuriousPower = 10.0*log10(maxSpuriousPower);
figure;
plot((1:PFBarms),maxSpuriousPower(1:PFBarms,1)); box on; grid on; axis([1 PFBarms -100 -50]); title('Max Spurious Power'); xlabel('Impulse Location'); ylabel('Max Spurious Power Component (dB)');

totalSpuriousPower = 10.0*log10(totalSpuriousPower);
figure;
plot((1:PFBarms),totalSpuriousPower(1:PFBarms,1)); box on; grid on; axis([1 PFBarms -100 -50]); title('Total Spurious Power'); xlabel('Impulse Location'); ylabel('Total Spurious Power (dB)');

averageSpuriousPower = 10.0*log10(averageSpuriousPower);
figure;
plot((1:blockLength),averageSpuriousPower(1:blockLength,1)); box on; grid on; axis([1 blockLength -150 -50]); title('Average Spurious Power'); xlabel('Time Sample'); ylabel('Average Spurious Power (dB)');

fprintf('\n\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
