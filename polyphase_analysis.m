function out=polyphase_analysis(in,filt,block,step)

% Polyphase analysis filterbank with cyclic shift of data into FFT 
% to remove spectrum rotation in output data
% John Bunton CSIRO 2003
% John Bunton CSIRO 2016

% in = input data
% filt = prototype lowpas filter (length should be multiple of step)
% block = length of fft (prefilter length = length(filt)/block
%           if not the 'filt' is padded with zeros to a multiple of block
% step = increment along 'in' between output blocks
%           step=block for critical sampling
% out = output data 2D one dimension time the other frequency

phases=ceil(length(filt)/block);
f=(1:phases*block)*0;
f(1:length(filt))=filt;

nblocks=floor( (length(in)-length(f))/step);
fl=length(f);

%block=block*2;     % Interleaved filterbank 
%phases=phases/2;   %produces critically sampled outputs as well as
                    %intermediate frequency outputs

for k=0:nblocks-1
    temp=f.*in(1+step*k:fl+step*k);
      
    %index for cyclic shift of data to FFT to eliminate spectrum rotation
    index = (step*k - floor(step*k/block)*block); 
    temp=circshift(temp',index)';
    
    temp2=(1:block)*0;
    for m=0:phases-1
        temp2=temp2+temp(1+block*m:block*(m+1));
    end
    out(k+1,1:block)= fft(temp2); %temp2;%
end
%plot(real(temp2))
end
