function y = Noise_vocoder_rcoef(x,Fmin,Fmax,Nband,Env_Cutoff,Fs,rcoef)
%This function creates envelope vocoded speech using noise carriers with a
%desired interaural correlation. Similar version used by Baltzell et al. 
%(2020) JASA: https://doi.org/10.1121/10.0000812. In this version, to
%prevent temporal distortions, zero-phase Butterworth filters are used to
%extract bands

%x: mono input (to yield stereo output)
%Fmin: minimum frequency
%Fmax: maximum frequency
%Nband: number of frequency bands used for synthesis
%Env_Cutoff: envelope cutoff frequency 
%Fs: sampling rate of input x
%rcoef: desired interaural correlation of noise carriers (0 to 1)

%ex. y = Noise_vocoder_rcoef(x,80,8000,32,300,Fs,0.5)

%created by Luke Baltzell 09/04/19, modified 04/19/21

[nch,dim] = min(size(x));
if dim == 1
    x = x';
end
if nch > 1
    error('should be mono input')
end

lev_input = rms(x); %rms of input

%initialize running average
env_accumulate_LE = zeros(size(x));
env_accumulate_RE = zeros(size(x));
 
%get log-spaced cutoff frequencies
Fctfs = logspace(log10(Fmin),log10(Fmax),Nband+1);
[sos_coefs,g_coefs] = get_butter_coefs(Fctfs,Fs,4);

%filter into bands
sig_filt = zeros(length(x),Nband);
for fb = 1:Nband
    sig_filt(:,fb) = filtfilt(sos_coefs{fb},g_coefs{fb},x);
end

%Generate noise carriers following Hartman & Cho (2013) JASA
noise1 = pinknoise(length(x),1,[20 Fs/3],Fs);
noise2 = pinknoise(length(x),1,[20 Fs/3],Fs);

alpha = sqrt((rcoef + 1)/2);
beta  = sqrt(1 - alpha^2);

noise_LE = alpha*noise1 + beta*noise2;
noise_RE = alpha*noise1 - beta*noise2;

for fb = 1:Nband
    noise_filt_LE(:,fb) = filtfilt(sos_coefs{fb},g_coefs{fb},noise_LE);
    noise_filt_RE(:,fb) = filtfilt(sos_coefs{fb},g_coefs{fb},noise_RE);
end

%Extract low pass filtered envelope
[Benv,Aenv]=butter(4,Env_Cutoff/(Fs/2));
filt_env=filtfilt(Benv,Aenv,abs(hilbert(sig_filt)));

for nband = 1:Nband
    %combine envelope and carrier
    env_LE(:,nband)=filt_env(:,nband).*cos(angle(hilbert(noise_filt_LE(:,nband))));
    env_RE(:,nband)=filt_env(:,nband).*cos(angle(hilbert(noise_filt_RE(:,nband))));
    
    %refilter band
    env_filt_LE(:,nband)=filtfilt(sos_coefs{nband},g_coefs{nband},env_LE(:,nband));
    env_filt_RE(:,nband)=filtfilt(sos_coefs{nband},g_coefs{nband},env_RE(:,nband));
    
    %add to running sum
    env_accumulate_LE=env_accumulate_LE+env_filt_LE(:,nband);   
    env_accumulate_RE=env_accumulate_RE+env_filt_RE(:,nband);
end

%apply a rise/fall window of 10 ms to avoid transients
env_accumulate_LE = rampdamp(env_accumulate_LE,0.01,Fs);
env_accumulate_RE = rampdamp(env_accumulate_RE,0.01,Fs);

%define scale factors
scale_LE = lev_input/rms(env_accumulate_LE);
scale_RE = lev_input/rms(env_accumulate_RE);

%scale each ear to match input level
y(:,1) = scale_LE*env_accumulate_LE;
y(:,2) = scale_RE*env_accumulate_RE;

end
