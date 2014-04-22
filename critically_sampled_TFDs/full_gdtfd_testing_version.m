%-------------------------------------------------------------------------------
% full_gdtfd_testing_version: TESTING version, *not* memory efficient
%
% Syntax: [tfd,g]=full_gdtfd_testing_version(x,kern_type,kern_params)
%
% Inputs: 
%     x,kern_type,kern_params - 
%
% Outputs: 
%     [tfd,g] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 15-04-2014
%
% last update: Time-stamp: <2014-04-22 10:49:44 (otoolej)>
%-------------------------------------------------------------------------------
function [tfd,g]=full_gdtfd_testing_version(x,kern_type,kern_params)
if(nargin<2 || isempty(kern_type)), kern_type='cw'; end
if(nargin<3 || isempty(kern_params)), kern_params=10; end

DBplot=0;
DBmem=0;
DBverbose=1;


%---------------------------------------------------------------------
% 1. convert real-valued signal to analytic signal of length 2N
%---------------------------------------------------------------------
[z,N2,N,Nh]=get_analytic_signal(x);


%---------------------------------------------------------------------
% 2. generate time--lag signal function (for positive-lag values only)
%---------------------------------------------------------------------
if(DBmem)
    s=whos; 
    fprintf('TESTING: start: mem=%g MB\n',disp_bytes(sum([s.bytes]))); 
end

K=zeros(N,Nh+1); 
m=0:(Nh);
for n=0:N-1
    inp=mod(n+m,N2);  inn=mod(n-m,N2); 
    K(n+1,m+1)=z(inp+1).*conj( z(inn+1) );
end

if(DBmem)
    s=whos; 
    fprintf('TESTING: K: mem=%g MB\n',disp_bytes(sum([s.bytes]))); 
end


%-------------------------------------------------------------------------
% 3. multiply kernel and signal function in the Doppler-lag domain
%-------------------------------------------------------------------------
g=gen_Doppler_lag_kern(kern_type,kern_params,N);

if(DBmem)
    s=whos; 
    fprintf('TESTING: g: mem=%g MB\n',disp_bytes(sum([s.bytes]))); 
end

R=ifft( fft(K(:,m+1)).*g(:,m+1) );
clear K;
if(DBmem)
    s=whos; 
    fprintf('TESTING: R: mem=%g MB\n',disp_bytes(sum([s.bytes]))); 
end


%-------------------------------------------------------------------------
% 4.  Expand R for positive and negative lag values
%-------------------------------------------------------------------------
Rfull=zeros(N);
Rfull(:,m+1)=R; 
mb=1:(Nh-1);
Rfull(:,N-mb+1)=conj( Rfull(:,mb+1) );
clear R;
if(DBmem)
    s=whos; 
    fprintf('TESTING: Rfull: mem=%g MB\n',disp_bytes(sum([s.bytes]))); 
end



%---------------------------------------------------------------------
% 5. back to the time--frequency domain
%---------------------------------------------------------------------
tfd=fft( Rfull.' ).';
clear Rfull;


if(DBverbose)
    fprintf('TESTING: max. imaginary component=%g\n', max(imag(tfd(:))) );
end
tfd=real(tfd)./N;


if(DBmem)
    s=whos; 
    fprintf('TESTING: tfd (end): mem=%s\n',disp_bytes(sum([s.bytes]))); 
end


if(DBplot)
    figure(1); clf; 
    vtfd(tfd,real(x(1:N)));
end

