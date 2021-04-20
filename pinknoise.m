function y = pinknoise(n,filt_flg,flims,fs)
%generate pink noise for a vector with n samples. If desired, apply
%high-pass and low-pass filter, requiring a sample rate
%n = number of samples
%filt_flg: 0 if no filter, 1 if filter
%flims: desired high-pass and low-pass cutoffs
%fs: sampling rate required for filtering. Must be 44000 Hz or greater

%created by Luke Baltzell 04/19/21

if nargin == 1
    filt_flg = 0;
end

x = randn(1,n); %generate white noise
f = [1:floor(n/2)+1]; %vector of frequencies

X = fft(x); %take FFT
Xpink = X(1:length(f))./sqrt(f); %1/f in power is 1/sqrt(f) in magnitude
if mod(length(x),2) == 0
    Xpink = cat(2,Xpink,fliplr(Xpink(3:end)));
else
    Xpink = cat(2,Xpink,fliplr(Xpink(2:end)));
end

y = ifft(Xpink,'symmetric');

if filt_flg == 1
    if fs < 44000
        error('sampling rate too low')
    end
    %apply high-pass filter at 20 Hz
    [b,a] = butter(2,flims(1)/(fs/2),'high');
    y = filter(b,a,y);
    
    %apply low-pass filter at 20 kHz
    [b,a] = butter(2,flims(2)/(fs/2),'low');
    y = filter(b,a,y);
end

end






