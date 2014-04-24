%-------------------------------------------------------------------------------
% dec_di_gdtfd: decimated Doppler-independent kernel TFD
%
% Syntax: tfd=dec_di_gdtfd(x,lag_win_params,time_dec,freq_dec,Nfreq)
%
% Inputs: 
%     x,lag_win_params,time_dec,freq_dec,Nfreq - 
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
% last update: Time-stamp: <2014-04-24 11:17:12 (otoolej)>
%-------------------------------------------------------------------------------
function tfd=dec_di_gdtfd(x,lag_win_params,time_dec,freq_dec,Nfreq)
if(nargin<2 || isempty(lag_win_params)), lag_win_params={101,'hamm'}; end
if(nargin<3 || isempty(time_dec)), time_dec=[]; end
if(nargin<4 || isempty(freq_dec)), freq_dec=1; end
if(nargin<5 || isempty(Nfreq)), Nfreq=[]; end

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
% 2. generate time--lag signal function (for positive-lag values only)
%---------------------------------------------------------------------
[g2,P,Ph_floor,Nfreq]=gen_lag_kern(lag_win_params,N,Nfreq);
Nh=ceil(Nfreq/2);
Ph=ceil(P/2);


% check decimation parameters and return as sequence:
if(length(freq_dec)>1)
    error('Frequency decimation parameter should be scalar');
end
[n_seq,U]=check_dec_params_seq(time_dec,N,'time',0);
[k_seq,J,freq_dec]=check_dec_params_seq(freq_dec,Nfreq,'frequency',1);
Jh=ceil(J/2);  J_extend=2*Jh+2;


if(DBmem), s=whos; fprintf('start: mem=%s\n',disp_bytes(sum([s.bytes]))); end
tfd=zeros(U,J_extend); 
if(J~=Nfreq)
    m_real=1:Jh+1; m_imag=(J_extend-Jh):(J_extend);
else
    m_real=1:Ph; m_imag=(J-Ph+1):J;
end
dispVars(length(m_real),length(m_imag),m_real(end),m_imag(1));


m=0:(Ph-1); mb=1:(Ph-1);
for in=1:U
    n=n_seq(in);
    inp=mod(n+m,N2);  inn=mod(n-m,N2); 
    
    if(J~=Nfreq)
        % need to fold signal function (for decimation in frequency-direction):        
        K_time_slice=zeros(1,Nfreq);
        K_time_slice(m+1)=g2(m+1).*z(inp+1).*conj( z(inn+1) );
        K_time_slice(Nfreq-mb+1)=conj(K_time_slice(mb+1));

        R_fold=fold_vector_half(K_time_slice,J,Jh,freq_dec);
    else
        % or if not then do as usual:
        K_time_slice(m+1)=g2(m+1).*z(inp+1).*conj( z(inn+1) );

        R_fold=K_time_slice;
    end
    
    tfd(in,m_real)=real( R_fold );
    tfd(in,m_imag)=imag( R_fold );    
end

if(DBmem), s=whos; fprintf('K: mem=%s\n',disp_bytes(sum([s.bytes]))); end


%-------------------------------------------------------------------------
% 3.  Expand R for positive and negative lag values and DFT back to 
%     time--frequency domain
%-------------------------------------------------------------------------
if(J~=Nfreq)
    m=0:Jh; mb=1:(Jh-1);
else
    m=0:(Ph-1); mb=1:(Ph-1);
end
    

for n=0:2:(U-2)
    R_even_half=complex(tfd(n+1,m_real),tfd(n+1,m_imag));
    R_odd_half =complex(tfd(n+2,m_real),tfd(n+2,m_imag));    
    
    R_tslice_even=zeros(1,J);  R_tslice_odd=zeros(1,J);
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_odd(m+1) =R_odd_half(m+1);
    R_tslice_even(J-mb+1)=conj( R_even_half(mb+1) );
    R_tslice_odd(J-mb+1) =conj( R_odd_half(mb+1) );
    
    tfd_time_slice=fft( R_tslice_even+j.*R_tslice_odd );

    tfd(n+1,1:J)=real(tfd_time_slice);
    tfd(n+2,1:J)=imag(tfd_time_slice);
end

% one extra FFT if U is odd:
if(rem(U,2))
    R_even_half=complex(tfd(U,m_real),tfd(U,m_imag));
    R_tslice_even=zeros(1,J);
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_even(J-mb+1)=conj( R_even_half(mb+1) );
    tfd_time_slice=fft( R_tslice_even );

    tfd(U,1:J)=real(tfd_time_slice);
end
tfd=tfd(:,1:J);

if(DBmem), s=whos; fprintf('end: mem=%s\n',disp_bytes(sum([s.bytes]))); end

scale_factor=1/Nfreq;
tfd=tfd.*scale_factor;

if(DBtime), dispVars( toc(time_start) ); end


%---------------------------------------------------------------------
% END; testing and plotting
%---------------------------------------------------------------------
if(DBtest)
    if(DBtime),  time_start=tic; end
    tfd_test=full_gdtfd_testing_version(x,'pwvd',{P,lag_win_params{2:end}});
    if(DBtime),  dispVars( toc(time_start) ); end

    if(Nfreq==N)
        dispEE(tfd_test(n_seq+1,k_seq+1),tfd);        
    else
        b=N/Nfreq; 
        if( b==floor(b) )
            dispEE(tfd_test(n_seq+1,(b.*k_seq)+1),tfd./b);
        end
    end

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
end


