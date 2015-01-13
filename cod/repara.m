%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alejandro Santa-Cruz                                           %
%   Reconstrucción de señales de MoCap por análisis de Fourier.    %
%   Instituto Tecnológico de Buenos Aires 2011-2012                %
%   fichero repara.m                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [REC,d]=repara(Igap,Dgap,T2,vars)


% Datos Iniciales 

    t=[Igap.t Dgap.t];
    x=[Igap.x Dgap.x];
    Ixgap=Igap.x;
    Dxgap=Dgap.x;
    NI=length(Ixgap);
    ND=length(Dxgap);
    N=length(T2);
    gaps=[Igap.t(end)+1:Dgap.t(1)-1];
    gapL=length(gaps);
    
    REC=zeros(1,N);
    REC(1:NI)=Igap.x;
    REC(N-ND+1:end)=Dgap.x;
    
    % contador de figuras
    f=1;   
    
    
%   CARACTERIZACION
%%%%%%%%%%%%%%%%%%%%%%%%%

% detectar frecuencias dominantes   
    pffdt= ND;                  pffit= NI;
    fftdx = fft(Dxgap,pffdt);   fftix = fft(Ixgap,pffit);
    nud = ceil((pffdt+1)/2);    nui = ceil((pffit+1)/2);
    fftdx = fftdx(1:nud);       fftix = fftix(1:nui);
    mxd = abs(fftdx)/ND;        mxi = abs(fftix)/NI; 
    mxd = mxd.^2;               mxi = mxi.^2; 
    if rem(pffdt, 2)                     % odd pfft excludes Nyquist point
      mxd(2:end) = mxd(2:end)*2;
    else
      mxd(2:end -1) = mxd(2:end -1)*2;
    end
    if rem(pffit, 2)                     % odd pfft excludes Nyquist point
      mxi(2:end) = mxi(2:end)*2;
    else
      mxi(2:end -1) = mxi(2:end -1)*2;
    end
    fd = (0:nud-1)*100/pffdt;           fi = (0:nui-1)*100/pffit; 
    % detectar picos maximos en ambos lados
    [maxi, mini] = peakdet(mxi, 0.05, fi);
    [maxd, mind] = peakdet(mxd, 0.05, fd);
    Ihay=0;Dhay=0;
    if length(maxi) > 0
        imaxind=find(maxi(:,2)>=maxi(1,2)*0.075);
        imaxind=find((100./maxi(imaxind)) < NI).';
        if not(isempty(imaxind))
            Ihay=1;
            ifreq=round(100./maxi(imaxind));
        end
    end
    if length(maxd) > 0
        dmaxind=find(maxd(:,2)>=maxd(1,2)*0.075);
        dmaxind=find((100./maxd(dmaxind))< ND).';
        if not(isempty(dmaxind))
            Dhay=1;
            dfreq=round(100./maxd(dmaxind));
        end
    end

    
    margen=30;
% Establecer un margen en los extremos
    if NI not(0);, iii=NI-margen;  if NI <= margen, iii=1;     end, end
    if ND not(0);, ddd=margen;     if ND <= margen, ddd=ND;    end, end
% Datos estadisticos de ambas partes
    iv=var(Ixgap(iii:end));
    im=median(abs(Ixgap(iii:end)-median(Ixgap(iii:end))));
    dv=var(Dxgap(1:ddd));
    dm=median(abs(Dxgap(1:ddd)-median(Dxgap(1:ddd))));
    difI=max(Ixgap(iii:end))-min(Ixgap(iii:end));
    difD=max(Dxgap(1:ddd))-min(Dxgap(1:ddd));


    
    
%   DECISION y RECONSTRUCCION
%%%%%%%%%%%%%%%%%%%%%%%%%
    % Condiciones y  Elementos de decision
    Freqs = (Ihay || Dhay);
    MayorE =   (difI >= vars(1) || difD >= vars(1)); %ind1    
    MenorE = not(difI < vars(2) || difD <  vars(2)); %ind2     
    MayorN =   (difI >= vars(3) || difD >= vars(3)); %ind3     
    MenorN = not(difI < vars(4) || difD <  vars(4)); %ind4
    CFou1 =    Freqs   && MayorE && MenorE;
    CFou2 = not(Freqs) && MayorN && MenorN;
    CondFourier = (CFou1 || CFou2);

    % CUB        
    if not(CondFourier)
        cub=interp1(t,x,T2,'*cubic');
        REC=cub;
        d='CUB';
        met=2;
    % FOU
    elseif CondFourier
        % Reordenar la señal: F(1)..F(290)   []   F(350)..F(800)
        %               ----> F(350)...F(800) F(1)...F(290)   []
        xr=[Dgap.x Igap.x];
        % Clacular F y S por edft
        W=ones(1,N);
        for it=1:2
            [F,S]=edft(xr,N,W);
            W=S;
        end
        % Calcular salida de EDFT
        sub=real(ifft(F));
        % Reordenar la señal
        fD=sub(1:ND);
        fI=sub(ND+1:ND+NI);
        nou=sub(ND+NI+1:end);
        % SUAVIZAR con fourier
        Y = fft(nou,gapL);
        fY=real(ifft(Y));
        sfou=fastsmooth(fY,5,2,1);
        % AJUSTAR desplazamiento
        off1=(Ixgap(end) - sfou(1));
        off2=(Dxgap(1) - sfou(end));
        fou=zeros(1,gapL);
        for i=1:gapL
            fou(i)=sfou(i) + off1*(gapL-i+1)/gapL + off2*(i-1)/gapL;       
        end
        
        REC(NI+1:N-ND)=fou;
        d='FOU';
        met=1;
    end