%-------------------------------------------------------------------------------
% di_gdtfd: TFD with Doppler-independent (DI) kernel
%           g[l,m]=g₂[m]
%
% Syntax: tfd=di_gdtfd(x,lag_win_params,Nfreq)
%
% Inputs: 
%     x,lag_win_params,Nfreq - 
%
% Outputs: 
%     tfd - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 16-04-2014
%
% last update: Time-stamp: <2014-04-22 11:53:26 (otoolej)>
%-------------------------------------------------------------------------------
function tfd=di_gdtfd(x,lag_win_params,Nfreq)
if(nargin<2 || isempty(lag_win_params)), lag_win_params={101,'hamm'}; end
if(nargin<3 || isempty(Nfreq)), Nfreq=[]; end


DBplot=1;
DBmem=1;
DBcompare=0;
DBtest=1;


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

if(DBmem), s=whos; fprintf('start: mem=%s\n',disp_bytes(sum([s.bytes]))); end
tfd=zeros(N,Nfreq); 
m_real=1:Ph; m_imag=(Nfreq-Ph+1):Nfreq;
dispVars(length(m_real),length(m_imag),m_imag(1));

m=0:(Ph-1);
for n=0:N-1
    inp=mod(n+m,N2);  inn=mod(n-m,N2); 
    K_time_slice=g2(m+1).*z(inp+1).*conj( z(inn+1) );
    
    tfd(n+1,m_real)=real( K_time_slice );
    tfd(n+1,m_imag)=imag( K_time_slice );    
end


if(DBcompare)
    K=zeros(N,Nh);
    for n=0:N-1
        inp=mod(n+m,N2);  inn=mod(n-m,N2); 
        K(n+1,m+1)=g2(m+1).*z(inp+1).*conj( z(inn+1) );
    end
    Ktest=complex(tfd(:,m_real),tfd(:,m_imag));
    dispEE(Ktest,K(:,m+1));
end
if(DBmem), s=whos; fprintf('K: mem=%s\n',disp_bytes(sum([s.bytes]))); end


%-------------------------------------------------------------------------
% 3.  Expand R for positive and negative lag values and DFT back to 
%     time--frequency domain
%-------------------------------------------------------------------------
m=0:(Ph-1); mb=1:(Ph-1);
for n=0:2:N-2
    R_even_half=complex(tfd(n+1,m_real),tfd(n+1,m_imag));
    R_odd_half =complex(tfd(n+2,m_real),tfd(n+2,m_imag));    
    
    R_tslice_even=zeros(1,Nfreq);  R_tslice_odd=zeros(1,Nfreq);
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_odd(m+1) =R_odd_half(m+1);
    R_tslice_even(Nfreq-mb+1)=conj( R_even_half(mb+1) );
    R_tslice_odd(Nfreq-mb+1) =conj( R_odd_half(mb+1) );
    
    tfd_time_slice=fft( R_tslice_even+j.*R_tslice_odd );

    tfd(n+1,:)=real(tfd_time_slice);
    tfd(n+2,:)=imag(tfd_time_slice);
end

if(DBmem), s=whos; fprintf('tfd: mem=%s\n',disp_bytes(sum([s.bytes]))); end


if(DBcompare)
    Rfull=zeros(N,Nfreq);
    Rfull(:,m+1)=K(:,m+1); 
    Rfull(:,Nfreq-mb+1)=conj( Rfull(:,mb+1) );

    tfd_test=fft( Rfull.' ).';
    dispEE(tfd_test,tfd);
end


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
        dispEE(tfd_test,tfd);
    else
        b=N/Nfreq; 
        if( b==floor(b) )
            dispEE(tfd_test(:,1:b:end),tfd./b);
        end
    end
end
if(DBplot)
    figure(1); clf; 
    vtfd(tfd,real(x(1:N)));
    
    figure(9); clf; hold all;
    subplot(211); hold all; plot(sum(tfd')'); plot( abs(z(1:N)).^2 );
    subplot(212); hold all; plot(sum(tfd')' - abs(z(1:N)).^2 );    
end

