%-------------------------------------------------------------------------------
% nonsep_gdtfd: Time-frequency distribution (quadratic class) with non-separable kernel
%
% Syntax: [tfd,g]=nonsep_gdtfd(x,kern_type,kern_params)
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
% Started: 14-04-2014
%
% last update: Time-stamp: <2014-05-02 18:54:48 (otoolej)>
%-------------------------------------------------------------------------------
function tfd=nonsep_gdtfd(x,kern_type,kern_params)
if(nargin<2 || isempty(kern_type)), kern_type='cw'; end
if(nargin<3 || isempty(kern_params)), kern_params=10; end


DBplot=0;
DBmem=1;
DBcompare=0;
DBtest=0;
DBtime=1;


if(DBtime), time_start=tic; end
%---------------------------------------------------------------------
% 1. convert real-valued signal to analytic signal of length 2N
%---------------------------------------------------------------------
[z,N2,N,Nh]=get_analytic_signal(x);


%---------------------------------------------------------------------
% 2. generate time--lag signal function (for positive-lag values only)
%---------------------------------------------------------------------
if(DBmem), s=whos; fprintf('start: mem=%s\n',disp_bytes(sum([s.bytes]))); end
tfd=zeros(N,N); 
m_real=1:Nh; m_imag=(Nh+1):N;

m=0:(Nh-1);
for n=0:N-1
    inp=mod(n+m,N2);  inn=mod(n-m,N2); 
    K_time_slice=z(inp+1).*conj( z(inn+1) );
    
    tfd(n+1,m_real)=real( K_time_slice );
    tfd(n+1,m_imag)=imag( K_time_slice );    
end

if(DBcompare)
    for n=0:N-1
        inp=mod(n+m,N2);  inn=mod(n-m,N2); 
        K(n+1,m+1)=z(inp+1).*conj( z(inn+1) );
    end
    Ktest=complex(tfd(:,m_real),tfd(:,m_imag));
    dispEE(Ktest,K);
end


if(DBmem), s=whos; fprintf('K: mem=%s\n',disp_bytes(sum([s.bytes]))); end


%-------------------------------------------------------------------------
% 3. multiply kernel and signal function in the Doppler-lag domain
%-------------------------------------------------------------------------
for m=0:Nh-1
    g_lag_slice=gen_Doppler_lag_kern(kern_type,kern_params,N,m+1);
    
    R_lag_slice=ifft( fft( complex(tfd(:,m_real(m+1)),tfd(:,m_imag(m+1))) ).*g_lag_slice );

    
    tfd(:,m_real(m+1))=real( R_lag_slice );
    tfd(:,m_imag(m+1))=imag( R_lag_slice );    
end

if(DBcompare)
    g=gen_Doppler_lag_kern(kern_type,kern_params,N);
    R=ifft( fft(K).*g(:,1:Nh) );
    Rtest=complex(tfd(:,m_real),tfd(:,m_imag));
    dispEE(Rtest,R);
    clear g;
end


if(DBmem), s=whos; fprintf('R: mem=%s\n',disp_bytes(sum([s.bytes]))); end




%-------------------------------------------------------------------------
% 4.  Expand R for positive and negative lag values and DFT back to 
%     time--frequency domain
%-------------------------------------------------------------------------
m=0:(Nh-1); mb=1:(Nh-1);
for n=0:2:N-2
    R_even_half=complex(tfd(n+1,m_real),tfd(n+1,m_imag));
    R_odd_half =complex(tfd(n+2,m_real),tfd(n+2,m_imag));    
    
    R_tslice_even=zeros(1,N);  R_tslice_odd=zeros(1,N);
    R_tslice_even(m+1)=R_even_half(m+1);
    R_tslice_odd(m+1) =R_odd_half(m+1);
    R_tslice_even(N-mb+1)=conj( R_even_half(mb+1) );
    R_tslice_odd(N-mb+1) =conj( R_odd_half(mb+1) );
    
    tfd_time_slice=fft( R_tslice_even+j.*R_tslice_odd );

    tfd(n+1,:)=real(tfd_time_slice);
    tfd(n+2,:)=imag(tfd_time_slice);
end

if(DBmem), s=whos; fprintf('tfd: mem=%s\n',disp_bytes(sum([s.bytes]))); end

    
if(DBcompare)
    Rfull=zeros(N);
    Rfull(:,m+1)=R; 
    Rfull(:,N-mb+1)=conj( Rfull(:,mb+1) );

    tfd_test=fft( Rfull.' ).';
    dispEE(tfd_test,tfd);
end

if(DBmem), s=whos; fprintf('end: mem=%s\n',disp_bytes(sum([s.bytes]))); end
tfd=tfd./N;


if(DBtime), dispVars( toc(time_start) ); end

if(DBtest)
    if(DBtime), time_start=tic; end
    tfd_test=full_gdtfd_testing_version(x,kern_type,kern_params);
    if(DBtime), dispVars( toc(time_start) ); end
    dispEE(tfd_test,tfd);
end
if(DBplot)
    figure(1); clf; 
    vtfd(tfd,real(x(1:N)));
end





