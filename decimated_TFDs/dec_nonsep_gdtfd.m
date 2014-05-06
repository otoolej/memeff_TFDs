%-------------------------------------------------------------------------------
% dec_nonsep_gdtfd: decimated TFD with non-separable kernel ρ[an,bk], where a,b are
% integer values
% 
%
% Syntax: tfd=dec_nonsep_gdtfd(x,kern_type,kern_params)
%
% Inputs: 
%     x,kern_type,kern_params - 
%
% Outputs: 
%     tfd - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 23-04-2014
%
% last update: Time-stamp: <2014-05-01 13:20:22 (otoolej)>
%-------------------------------------------------------------------------------
function tfd=dec_nonsep_gdtfd(x,kern_type,kern_params,time_dec,freq_dec)
if(nargin<2 || isempty(kern_type)),   kern_type='cw'; end
if(nargin<3 || isempty(kern_params)), kern_params=10; end
if(nargin<4 || isempty(time_dec)),    time_dec=1; end
if(nargin<5 || isempty(freq_dec)),    freq_dec=1; end


DBplot=1;
DBmem=1;
DBtest=1;
DBtime=0;
DBverbose=1;


if(DBtime), time_start=tic; end
%---------------------------------------------------------------------
% 1. convert real-valued signal to analytic signal of length 2N
%---------------------------------------------------------------------
[z,N2,N,Nh]=get_analytic_signal(x);


%---------------------------------------------------------------------
% 2. check decimation parameters 
%---------------------------------------------------------------------
if(length(time_dec)>1 || length(freq_dec)>1)
    error('Frequency and time decimation parameters should be scalar');
end
[n_seq,L,time_dec]=check_dec_params_seq(time_dec,N,'time',1);
[k_seq,J,freq_dec]=check_dec_params_seq(freq_dec,N,'frequency',1);
Jh=ceil(J/2);  J_extend=2*Jh+2;

if(DBmem), s=whos; fprintf('start: mem=%s\n',disp_bytes(sum([s.bytes]))); end

tfd=zeros(L,J_extend); 
m_real=1:Jh+1; m_imag=(J_extend-Jh):(J_extend);

if(DBmem), s=whos; fprintf('declare TFD: mem=%s\n',disp_bytes(sum([s.bytes]))); end


%---------------------------------------------------------------------
% 3. generate time--lag signal function (for positive-lag values only)
%---------------------------------------------------------------------
n=0:N-1; mb=1:Jh-1;
for m=0:Jh

    %------------------------------------
    % Fold in the lag direction
    %------------------------------------
    af_lag_slice=zeros(N,1);    
    for p=0:freq_dec-1
        mmod=(p*J+m);
        g_lag_slice=gen_Doppler_lag_kern(kern_type,kern_params,N,mmod+1);

        if(mmod<=Nh)
            inp=mod(n+mmod,N2); inn=mod(n-mmod,N2);
            K_lag_slice=z(inp+1).*conj(z(inn+1));
        else
            inp=mod(n+N-mmod,N2); inn=mod(n-N+mmod,N2);
            K_lag_slice=conj(z(inp+1)).*z(inn+1);            
        end
        
        af_lag_slice=af_lag_slice + fft(K_lag_slice).*g_lag_slice;
    end

    %------------------------------------
    % Fold in the Doppler direction
    %------------------------------------
    af_lag_fold=fold_vector_full(af_lag_slice,L,time_dec);
    
    %------------------------------------
    % DFT the lag slice to the time--lag 
    % domain.
    %------------------------------------
    R_lag_slice=ifft(af_lag_fold);
    tfd(:,m_real(m+1))=real(R_lag_slice);
    tfd(:,m_imag(m+1))=imag(R_lag_slice);    
end


if(DBmem), s=whos; fprintf('R: mem=%s\n',disp_bytes(sum([s.bytes]))); end

%-------------------------------------------------------------------------
% 4.  Expand R for positive and negative lag values and DFT back to 
%     time--frequency domain
%-------------------------------------------------------------------------
m=0:Jh; mb=1:(Jh-1);
for n=0:2:(L-2)
    R_even_half=complex(tfd(n+1,m_real),tfd(n+1,m_imag));
    R_odd_half =complex(tfd(n+2,m_real),tfd(n+2,m_imag));
    
    R_tslice_even=zeros(J,1);  R_tslice_odd=zeros(J,1);
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_odd(m+1) =R_odd_half(m+1);
    R_tslice_even(J-mb+1)=conj( R_even_half(mb+1) );
    R_tslice_odd(J-mb+1) =conj( R_odd_half(mb+1) );
    
    tfd_time_slice=fft( R_tslice_even+j.*R_tslice_odd );

    tfd(n+1,1:J)=real(tfd_time_slice);
    tfd(n+2,1:J)=imag(tfd_time_slice);
end


% one extra FFT if L is odd:
if(rem(L,2))
    R_even_half=complex(tfd(L,m_real),tfd(L,m_imag));
    
    R_tslice_even=zeros(J,1);  
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_even(J-mb+1)=conj( R_even_half(mb+1) );
    
    tfd_time_slice=fft( R_tslice_even+j.*R_tslice_odd );

    tfd(L,1:J)=real(tfd_time_slice);
end

tfd=tfd(:,1:J);

if(DBmem), s=whos; fprintf('end: mem=%s\n',disp_bytes(sum([s.bytes]))); end

scale_factor=1/N;
tfd=tfd.*scale_factor;


if(DBtime), dispVars( toc(time_start) ); end

%---------------------------------------------------------------------
% END; testing and plotting
%---------------------------------------------------------------------
if(DBtest)
    if(DBtime),  time_start=tic; end
    tfd_test=full_gdtfd_testing_version(x,kern_type,kern_params);
    if(DBtime),  dispVars( toc(time_start) ); end

    dispEE(tfd_test(n_seq+1,k_seq+1),tfd);
end
if(DBverbose)
    fprintf('size=%dx%d; max=%g; total energ:%d\n', size(tfd,1), size(tfd,2), max(tfd(:)), ...
            sum(tfd(:)));
end
if(DBplot)
    figure(1); clf; 
    vtfd(tfd,real(x(1:N)));
    
    figure(9); clf; hold all;
    subplot(211); hold all; plot(sum(tfd')'); plot( abs(z(1:N)).^2 );
    subplot(212); hold all; plot(sum(tfd')' - abs(z(n_seq+1)).^2 );    
    
    figure(10); clf; hold all;
    Z=fft(z);
    subplot(211); hold all; plot(sum(tfd)); plot( abs(Z(1:N)).^2./(2*N) );
    subplot(212); hold all; plot(sum(tfd)' - abs(Z(k_seq+1)).^2./(2*N) );    
end
