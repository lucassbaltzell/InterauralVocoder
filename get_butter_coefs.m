function [sos,g] = get_butter_coefs(cutoffs,Fs,N)
%rather than the transfer function [b,a], we obtain the second-order
%section representation for each band, which is more numerically stable for
%higher-order filters

%created by Luke Baltzell 04/19/21

if nargin == 2
    N = 8; %48 dB/octave
end

Nq = Fs/2;
sos = cell(length(cutoffs)-1,1);
g = cell(length(cutoffs)-1,1);
for i = 1:length(cutoffs)-1
    Wn = [cutoffs(i) cutoffs(i+1)]/Nq;
    [z,p,k] = butter(N,Wn);
    [sos{i}, g{i}] = zp2sos(z,p,k);
end