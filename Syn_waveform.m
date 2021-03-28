% Create synthetic waveforms
% Perform Time-Frequency Analysis (Wigner distribution)

rdist    = [100 150 200 250];   %  Receiver Distances
fmax     = 5.00;                %  Maximum freq
fmin     = 0.05;                %  Minimum freq
vmin     = 1.9;                 %  Minimum Phase Velocity
vmax     = 3.5;                 %  Maximum phase velocity

freq     = linspace(0,10,1024);
df       = freq(2)-freq(1);
imin     = fix(fmin/df);
imax     = fix(fmax/df);
frqsel   = (imin:imax)*df;
nfsel    = length(frqsel);
dv       = (vmax-vmin)/(imax-imin);
pvelsel  = vmax-(imin:imax)*dv;
nf       = length(freq);
deltaf   = max(freq)/nf;

dt       = 1/(2*max(freq));
tmax     = 2*nf*dt;
nt       = fix(tmax/dt);
nt       = fix(nt/2);
nt       = nt*2;
nsample  = nt;
nr       = 4;
tm       = (0:nt-1)*dt;

im       = sqrt(-1);
nfw      = fix(nfsel/20);
u1_fft   = zeros(nt,4);
wvfr     = zeros(nt,5);
iww      = 0;
for iw=2:nf
   if ( freq(iw)>frqsel(1) & freq(iw)<frqsel(nfsel) )       
      iww=iww+1;
      if iww<=nfw
         rw         = 0.5*(1-cos(pi*(iww-1)/nfw));
      elseif iww>=(nfsel-nfw)
         rw         = 0.5*(1-cos(pi*(nfsel-iww-1)/nfw));
      else
         rw = 1;
      end
      %rw=1;
      arg = -im * 2 * pi * frqsel(iww) ;
      u1_fft(iw,1:nr) = rw*exp( arg * ( rdist(1:nr) / pvelsel(iww) ));
   else
      u1_fft(iw,1:nr) = 0.0;
   end
end
phs=(atan2(imag(u1_fft),real(u1_fft)));
amp=abs(u1_fft);
u1_fft = amp.* (cos(phs) + im * sin(phs) );
for ir=1:4
for iw=1:nf-1
   u1_fft(nt-iw+1,ir) = conj(u1_fft(iw+1,ir));
end
   u1_fft(1,ir) = 0.0 ;
   u1_fft(nt/2,ir) = 0.0 ;
end
ibeg=1;
iend=nt;
pltmod='a';
scltrc=1;
xlabel('Time (Sec)','Fontsize',8);      
wvfr(:,1) = real(ifft(u1_fft(:,1)));
nsample      = nt;
ibeg         = 1; 
iend         = nsample; 
subplot(2,1,1), plot(tm,wvfr(:,1));
xlabel('Time (sec)')
ylabel('Receiver No')
axis([min(tm) max(tm) min(wvfr(:,1)) max(wvfr(:,1)) ]) 
[y,f,t,p] = spectrogram(wvfr(:,1),128,16,128);
subplot(2,1,2), imagesc(t,f,10*log10(abs(p)));
xlabel('Time (Sec)');      
ylabel('Frequency (Hz)');      
