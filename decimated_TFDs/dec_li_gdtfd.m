%-------------------------------------------------------------------------------
% dec_li_gdtfd: decimated lag-independent kernel TFD ρ[an,kᵢ], where a is an integer
% value and kᵢ is the set kᵢ ={ kᵢ | 1≤i≤V }, with 0≤kᵢ≤N-1
%
% Syntax: tfd=dec_li_gdtfd(x,dopp_win_params,time_dec,freq_dec,Ntime)
%
% Inputs: 
%     x,dopp_win_params,time_dec,freq_dec,Ntime - 
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
% last update: Time-stamp: <2014-05-01 14:35:22 (otoolej)>
%-------------------------------------------------------------------------------
function tfd=dec_li_gdtfd(x,dopp_win_params,time_dec,freq_dec,Ntime)
if(nargin<2 || isempty(dopp_win_params)), dopp_win_params={51,'hann'}; end
if(nargin<3 || isempty(time_dec)), time_dec=1; end
if(nargin<4 || isempty(freq_dec)), freq_dec=1; end
if(nargin<5 || isempty(Ntime)), Ntime=[]; end


DBplot=1;
DBmem=1;
DBtest=1;
DBtime=1;
DBverbose=1;


if(DBtime), time_start=tic; end
%---------------------------------------------------------------------
% 1. convert real-valued signal to analytic signal of length 2N
%---------------------------------------------------------------------
[z,N2,N,Nh]=get_analytic_signal(x);


%---------------------------------------------------------------------
% 2. generate time--lag signal function (for positive-lag values only)
%---------------------------------------------------------------------
[G1,Q,Ntime]=gen_Doppler_kern(dopp_win_params,N,Ntime);
Qh=ceil(Q/2);


% check decimation parameters and return as sequence:
if(length(time_dec)>1)
    error('Frequency decimation parameter should be scalar');
end
[n_seq,L,time_dec]=check_dec_params_seq(time_dec,Ntime,'time',1);
[k_seq,V]=check_dec_params_seq(freq_dec,N,'frequency',0);
Lh=ceil(L/2);  L_extend=2*Lh+2;


if(DBmem), s=whos; fprintf('start: mem=%s\n',disp_bytes(sum([s.bytes]))); end
tfd=zeros(L_extend,V); 
if(L~=Ntime)
    l_real=1:(Lh+1); l_imag=(L_extend-Lh):L_extend;
else
    l_real=1:Qh; l_imag=(L-Qh+1):L;
end
dispVars(length(l_real),length(l_imag),l_real(end),l_imag(1));


Z=fft(z);
l=0:(Qh-1); lb=1:(Qh-1);
for ik=1:V
    k=k_seq(ik);
    inp=mod(k+l,N2);  inn=mod(k-l,N2); 
    inpN=mod(k+l+N,N2);  innN=mod(k-l-N,N2);     
    
    if(L~=Ntime)
        K_doppler_slice=zeros(1,Ntime);
        K_doppler_slice(l+1)=G1(l+1).*( Z(inp+1).*conj( Z(inn+1) )  + ...
                                        Z(inpN+1).*conj( Z(innN+1) ) );
        K_doppler_slice(Ntime-lb+1)=conj( K_doppler_slice(lb+1) );
        
        R_fold=fold_vector_half(K_doppler_slice./2,L,Lh,time_dec);
    else
        
        K_doppler_slice=G1(l+1).*( Z(inp+1).*conj( Z(inn+1) )  + ...
                                   Z(inpN+1).*conj( Z(innN+1) ) );
        
        R_fold=K_doppler_slice;
    end
    
    tfd(l_real,ik)=real( R_fold );
    tfd(l_imag,ik)=imag( R_fold );    
end


if(DBmem), s=whos; fprintf('K: mem=%s\n',disp_bytes(sum([s.bytes]))); end



%-------------------------------------------------------------------------
% 3.  Expand R for positive and negative lag values and DFT back to 
%     time--frequency domain
%-------------------------------------------------------------------------
if(L~=Ntime)
    l=0:Lh; lb=1:(Lh-1);
else
    l=0:(Qh-1); lb=1:(Qh-1);
end

for k=0:2:(V-2)
    R_even_half=complex(tfd(l_real,k+1),tfd(l_imag,k+1));
    R_odd_half =complex(tfd(l_real,k+2),tfd(l_imag,k+2));    
    
    R_tslice_even=zeros(L,1);  R_tslice_odd=zeros(L,1);
    R_tslice_even(l+1)=R_even_half(l+1);
    R_tslice_odd(l+1) =R_odd_half(l+1);
    R_tslice_even(L-lb+1)=conj( R_even_half(lb+1) );
    R_tslice_odd(L-lb+1) =conj( R_odd_half(lb+1) );
    
    tfd_freq_slice=ifft( R_tslice_even+j.*R_tslice_odd );

    tfd(1:L,k+1)=real(tfd_freq_slice);
    tfd(1:L,k+2)=imag(tfd_freq_slice);
end

% one extra FFT if V is odd
if(rem(V,2))
    R_even_half=complex(tfd(l_real,V),tfd(l_imag,V));
    
    R_tslice_even=zeros(L,1);  
    R_tslice_even(l+1)=R_even_half(l+1);
    R_tslice_even(L-lb+1)=conj( R_even_half(lb+1) );
    
    tfd_freq_slice=ifft( R_tslice_even );

    tfd(1:L,V)=real(tfd_freq_slice);
end
tfd=tfd(1:L,:);

if(DBmem), s=whos; fprintf('tfd (end): mem=%s\n',disp_bytes(sum([s.bytes]))); end


if(time_dec==1)
    scale_factor=1/(2*N);
else
    scale_factor=1/(N*time_dec);
end
tfd=tfd.*scale_factor;

if(DBtime), dispVars( ['time: ' num2str(toc(time_start))] ); end


%---------------------------------------------------------------------
% END; testing and plotting
%---------------------------------------------------------------------
if(DBtest)
    if(DBtime), time_start=tic; end
    tfd_test=full_gdtfd_testing_version(x,'swvd',{Q,dopp_win_params{2:end}});
    if(DBtime), dispVars( ['testing time: ' num2str(toc(time_start))] ); end
    if(Ntime==N)
        dispEE(tfd_test(n_seq+1,k_seq+1),tfd);
    else
        a=N/Ntime; 
        if( a==floor(a) )
            dispEE(tfd_test((a.*n_seq)+1,k_seq+1),tfd./a);
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
    Z=fft(z);
    subplot(211); hold all; plot(sum(tfd)); plot( abs(Z(1:N)).^2./(2*N) );
    subplot(212); hold all; plot(sum(tfd)' - abs(Z(k_seq+1)).^2./(2*N) );    
end
