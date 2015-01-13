%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alejandro Santa-Cruz                                           %
%   Reconstrucción de señales de MoCap por análisis de Fourier.    %
%   Instituto Tecnológico de Buenos Aires 2011-2012                %
%   fichero procesa.m                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [data,sdif,sdcs,o,gapL]=procesa(to,xo,v)
%
%   INPUT
%   xo(to): señal uniforme completa
%   v: modo verbose
%   OUTPUT
%   sdif: suma de las diferencias aboslutas (MAE)
%   sdcs: suma de los cuadrados de las diferencias (RMSE)
%   data: datos de la caracterizacion (Amplitud, varianza, desviacion de mediana)
%   o: observacion de la caracterizacion
%   gapL: Longitud del hueco
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,sdif,sdcs,o,gapL]=procesa(to,xo,v)
    if nargin == 0
        % algunos datos por defecto para hacer pruebas
        addpath('data');
        load xmaxi10.pts
        z = xmaxi10';
        marker=randi(34,1,1).';
        to = z(1,:);K=to(end);
        xo = z(1+marker,:);
        v=1;    % sacar resultados
    end
    
    
    %%% INTRO
    %%%%%%%%%%%%%%%%%%%
    % contador de figuras
    f=1;
    K=to(end);
    
    gapI=1+randi(K-3,1,1);  % gapI: Inicio del gap
    gapL=50+randi(250,1,1);    % gapL: Longitud del gap >50
    gapF=gapI+gapL-1;       % gapF: Fin del gap
    if gapF>=K, gapF=K-2; end   % ajustar los extremos
    gaps = [gapI:gapF];
    gapL = length(gaps);
    % para graficar:
    xgap=xo(gaps);tgap=to(gaps);
    xnan=xo;xnan(gaps)=NaN;
    
%
% rangos observacion
% Parámetros: bs-boundaries; gapA-extremo1; gapZ-extremo2	
    bs=300;
    gapA = gaps(1) - bs;
    gapZ = gaps(end) + bs;
    if gapA < 1     % Comportamiento en los extremos
        gapA=1;
    elseif gapZ > K;
        gapZ=K;
    end
%                                     xrI    +     xrG     +   xrD     = N
% Datos, Rango para trabajar: xb = [extremo1 + gap1...gapN + extremo2]
    rangbnd=horzcat((gapA:gapI-1),(gapI:gapF),(gapF+1:gapZ));
    tb=to(rangbnd); 
    xb=xo(rangbnd);
    N=length(xb);
    fn=[-ceil((N-1)/2):floor((N-1)/2)]/N;
    
% Generar x non-uniform - Rango xb con gaps
    rangaps=horzcat((gapA:gapI-1),(gapF+1:gapZ));
    t=to(rangaps);
    x=xo(rangaps);
    Ixgap=xo(gapA:gapI-1);  xrI=length(Ixgap);
    Itgap=to(gapA:gapI-1);
    Dxgap=xo(gapF+1:gapZ);  xrD=length(Dxgap);
    Dtgap=to(gapF+1:gapZ);
    ND=length(Dxgap);           NI=length(Ixgap);    
%
% rangos precision
% Establecer un margen en los extremos
    margen=30;
    if NI not(0);, wi=NI-margen;  if NI <= margen, wi=1;     end, end
    if ND not(0);, wd=margen;     if ND <= margen, wd=ND;    end, end


% RECONSTRUCCIÓN INT SPLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spl=interp1(t,x,tgap,'*spline');    

% RECONSTRUCCIÓN INT CUBICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cub=interp1(t,x,tgap,'*cubic');

% RECONSTRUCCIÓN DE FOURIER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reordenar la señal: de F(1)..F(299) [] F(301)..F(800) ----> F(301)...F(800)F(1)...F(299) []
    rangdesord=horzcat((gapF+1:gapZ),(gapA:gapI-1));
    xr=xo(rangdesord);
% Parámetros: I-iteraciones;	
    I=2;
    xrD=length(Dxgap);
    xrI=length(Ixgap);
    % Clacular F y S por edft
    W=ones(1,N);
    for it=1:I
        [F,S]=edft(xr,N,W);
        W=S;
    end
    % Calcular salida de EDFT
    % sub1=10*log10(fftshift(abs(F).^2/N)); %PSD
    % sub2=20*log10(fftshift(abs(S)));  %PS
    % sub3=real(fftshift(F./S)/K);  %PSD/PS
    sub4=real(ifft(F));
% % Reordenar la señal
    fD=sub4(1:xrD);
    fI=sub4(xrD+1:xrD+xrI);
    nou=sub4(xrD+xrI+1:end);
    % SUAVIZAR con fourier
    Y = fft(nou,gapL);
    fY=real(ifft(Y));
    if gapL > 5
        sfou=fastsmooth(fY,5,2,1);
    else
        sfou=fY;
    end
    % AJUSTAR desplazamiento
    off1=(Ixgap(end) - sfou(1));
    off2=(Dxgap(1) - sfou(end));
    fou=zeros(1,gapL);
    for i=1:gapL
        fou(i)=sfou(i) + off1*(gapL-i+1)/gapL + off2*(i-1)/gapL;       
    end
    
    
    
    
%   OBSERVACIÓN  %
%%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimación de la calidad de reconstrucción
%    % valores para global        % errores de recons locales
    sdif.spl=sum(abs(xgap-spl)); lsae=(sdif.spl)./gapL;
    sdif.cub=sum(abs(xgap-cub)); lcae=(sdif.cub)./gapL;
    sdif.fou=sum(abs(xgap-fou)); lfae=(sdif.fou)./gapL;
    sdcs.spl=sum((xgap-spl).^2); lsms=sqrt((sdcs.spl)./gapL);
    sdcs.cub=sum((xgap-cub).^2); lcms=sqrt((sdcs.cub)./gapL);
    sdcs.fou=sum((xgap-fou).^2); lfms=sqrt((sdcs.fou)./gapL);
% * Método X mejor que método Y: Si MAE_X / MAE_Y < 0.8
% * Métodos equivalentes: Si 0.8 < MAE_X / MAE_Y < 1.2
% * Método Y mejor que método X:
    OKFou = lfae/lcae <= 0.8;
    OKCub = lfae/lcae >= 1.2;
    % equiv = lfae/lcae > 0.8  &&  lfae/lcae < 1.2;
    OKSpl = gapL <= 50;
    
    if OKSpl
        o.num=3;
        o.str='SPL';
    elseif not(OKFou) && not(OKCub)
        o.num=0;
        o.str='EQU';
    elseif OKFou
        o.num=1;
        o.str='FOU';
    elseif OKCub
        o.num=2;
        o.str='CUB';
    end

    
    
    
    
% CARACTERIZACIÓN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caracterzacion global

% FFTs anterior y posterior    
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
    % detectar picos maximos
    delta=0.05;
    [maxi, mini] = peakdet(mxi, delta, fi);
    [maxd, mind] = peakdet(mxd, delta, fd);
    Ihay=0;Dhay=0;
    if length(maxi) > 0
        imaxind=find(maxi(:,2)>=maxi(1,2)*0.075);
        imaxind=find((100./maxi(imaxind)) < bs).';
        if not(isempty(imaxind))
            Ihay=1;
            ifreq=round(100./maxi(imaxind));
        end
    end
    if length(maxd) > 0
        dmaxind=find(maxd(:,2)>=maxd(1,2)*0.075);
        dmaxind=find((100./maxd(dmaxind))< bs).';
        if not(isempty(dmaxind))
            Dhay=1;
            dfreq=round(100./maxd(dmaxind));
        end
    end

    
    % Caracterzacion local
    % Datos estadisticos de ambas partes
    data.iv=var(Ixgap(wi:end));
    data.im=median(abs(Ixgap(wi:end)-median(Ixgap(wi:end))));
    data.dv=var(Dxgap(1:wd));
    data.dm=median(abs(Dxgap(1:wd)-median(Dxgap(1:wd))));
    difI=max(Ixgap(wi:end))-min(Ixgap(wi:end));
    difD=max(Dxgap(1:wd))-min(Dxgap(1:wd));
    data.idst=difI;
    data.ddst=difD;
    data.iay=Ihay;
    data.day=Dhay;

    
%   ELEMENTOS DECISION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solamente INFORMATIVO y sacar graficos.
% Cambiar parametros de configuracion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dec=xo;

if gapL<=50
    dec(gaps)=spl;
    d.str='SPL';
    d.num=3;
else

% Condiciones Cubic y Fourier
    vars=[466,7,325,19];
    Freqs = (Ihay || Dhay);
    MayorE =   (difI >= vars(1) || difD >= vars(1)); %ind1
    menorE = not(difI < vars(2) || difD <  vars(2)); %ind2
    MayorN =   (difI >= vars(3) || difD >= vars(3)); %ind3
    menorN = not(difI < vars(4) || difD <  vars(4)); %ind4
    CF1 =    Freqs   && MayorE && menorE;
    CF2 = not(Freqs) && MayorN && menorN;
    CondFourier = (CF1 || CF2);
    
    if CondFourier
        dec(gaps)=fou;
        d.str='FOU';
        d.num=1;
    elseif not(CondFourier)
        dec(gaps)=cub;
        d.str='CUB';
        d.num=2;
    end 
end



%   
%   Algunas pruebas realizadas
%%%%%%%%%%%%%%%%%%%%%%

% Comparación de transformadas FFT vs EDFT
%
    nuds = ceil((N+1)/2);   fftdx = fft(xb,N);
    F2 = F(1:nuds);         fftdx = fftdx(1:nuds);
    edftfx = abs(F2)/N;     fx = abs(fftdx)/ND;
    edftfx = edftfx.^2;     fx = fx.^2;
    if rem(N, 2)
      edftfx(2:end) = edftfx(2:end)*2;
      fx(2:end) = fx(2:end)*2;
    else
      edftfx(2:end -1) = edftfx(2:end -1)*2;
      fx(2:end -1) = fx(2:end -1)*2;
    end
    ft = (0:nuds-1)*100/N;    


fprintf(['MoCap: [', num2str(K), '] Obs(i)=', o.str, ',  Hueco(',int2str(gapI),':',int2str(gapF), ')',...
    '\tTam=',int2str(gapL), '\t\tFreqs=',int2str((Ihay || Dhay)), '  Meds=',int2str(difI), ',' ,int2str(difD),'\n']);

    
    
    
if v==1
   
   
    dir = './TEST';
    if not(isdir(dir)),mkdir(dir);end
   
%
% DISPLAY Y GRAFICAS
%%%%%%%%%%%%%%%%%%%%

    % valores para graficas
    Xmin=min(xo);
    Xmax=max(xo);

% %
% % Comparación de la FFT con la EDFT
% %
    % figure(f);f=f+1;
    % set(gcf,'Visible','off');
    % clf
    % subplot(211)
    % plot(tb,xb,'k',tgap,xgap,'r')
    % % axis([gapA gapZ Xmin-50 Xmax+50])
    % ylabel(['Posicion'])
    % title(['Señal original y con un hueco'])
    % xlabel('Frames')
    % set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    % legend('Señal','Hueco','Location','NorthEast');
    % subplot(223)
    % semilogx(ft,fx,'b')
    % ylabel(['Unidad Arbitraria'])
    % title([' Power Spectral Density FFT'])
    % xlabel('Frecuencias')
    % set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    % legend('FFT','Location','NorthEast');
    % subplot(224)
    % semilogx(ft,fx,'b',ft,edftfx,'g')
    % ylabel(['Unidad Arbitraria'])
    % title([' Densidad de Potencia Espectral EDFT'])
    % xlabel('Frecuencias')
    % set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    % legend('FFT','EDFT','Location','NorthEast');
    % filenam = [dir,'/edft_x','_gaps(',num2str(gapI),'_',num2str(gapF),')'];
    % print(gcf,'-dpng','-r125',filenam);

% %
% % gráfica con toda la informacion de la reconstruccion y estudio
% %
    figure(f);f=f+1;
    set(gcf,'Visible','off');
    clf
    subplot(421)
    semilogx(fi,mxi)
    if Ihay
        hold on
        semilogx(maxi(imaxind,1), maxi(imaxind,2), 'r*')
        xlabel(['   Posible freq cada ',num2str(ifreq),' frames']);
    end
    title('    Detectar Frecuencias Dominantes')
    subplot(422)
    semilogx(fd,mxd)
    if Dhay
        hold on
        semilogx(maxd(dmaxind,1), maxd(dmaxind,2), 'r*')
        xlabel(['   Posible freq cada ',num2str(dfreq),' frames']);
    end
    title('    Detectar Frecuencias Dominantes')
    subplot(423)
    plot(Itgap(wi:end),Ixgap(wi:end),'k')
    % xlabel(['d: ',num2str(difI),',   v: ',num2str(data.iv),',   m: ',num2str(data.im)]);
    xlabel(['amp: ',num2str(difI)]);
    set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    title('Mediciones:')
    subplot(424)
    plot(Dtgap(1:wd),Dxgap(1:wd),'k')
    % xlabel(['d: ',num2str(difD),',   v: ',num2str(data.dv),',   m: ',num2str(data.dm)]);
    xlabel(['amp: ',num2str(difD)]);
    set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    title('Mediciones:')
    subplot(4,2,[5:6])
    plot(to,xnan,'k',tgap,xgap,'r')
    hold on
    plot(tgap,fou,'g',tgap,cub,'b',tgap,spl,'c')
    hold on
    plot(tb(1),xb(1),'b*',tb(end),xb(end),'b*')
    axis([gapA gapZ Xmin-50 Xmax+50])
    ylabel(['Posicion'])
    xlabel([' (', int2str(gapI),':',int2str(gapF),') T=',int2str(gapL)])
    set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    legend('Señal','Hueco','FOU','CUB','SPL','Location','BestOutside');
    subplot(4,2,[7:8])
    plot(to,xnan,'k',tgap,xgap,'r')
    hold on
    plot(to,dec,'k')
    hold on
    plot(tb(1),xb(1),'b*',tb(end),xb(end),'b*')
    ylabel(['Posicion'])
    xlabel(['Obs(i)= ', o.str, '    MAE_{FOU}= ',num2str(lfae), ',   MAE_{CUB}= ',num2str(lcae), ',   MAE_{SPL}= ',num2str(lsae)])
    set(gca,'xminortick','on','fontsize',11,'linewidth',1,'ticklength',[0.015 0.015])
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 20 15])
    filenam = [dir,'/edft_x','_gaps(',num2str(gapI),'_',num2str(gapF),')'];
    print(gcf,'-dpng','-r125',filenam);
    
end
    
    
    
    