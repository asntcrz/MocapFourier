%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alejandro Santa-Cruz                                           %
%   Reconstrucción de señales de MoCap por análisis de Fourier.    %
%   Instituto Tecnológico de Buenos Aires 2011-2012                %
%   fichero mocap_reconstr.m                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fichero principal de reconstruccion de los datos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tr,xr,resumen]=mocap_reconstr(t,x,p)
%
%   INPUT
%   t: valores del tiempo (pueden ser frames, por ej: tasa=100fps)
%   x: valores de la posicion del marcador
%   p: [OPCIONAL] parametros para la reconstruccion [p1,p2,p3,p4]
%   OUTPUT
%   tr: secuencia de tiempo completa, 1..t(end)
%   xr: secuencia de datos reconstruida
%   resumen: estructura con info acerca de las reconstrucciones
%       resumen.cs: #recs por metodo SPL
%       resumen.cc: #recs por metodo CUB
%       resumen.cf: #recs por metodo FOU
%       resumen.s:  tiempo en reconstruir el marcador entero
%
% Comprobacion básica de entradas
    if nargin == 2
        % Valores por defecto
        p=[466,7,325,19];
    elseif nargin == 3
        if size(p) ~= 4
            error('Error en los parametros','Debe haber 4 parametros');
        end
    else
        error('Error en la llamada','Mal uso de la llamada');
    end
    
    tic;
    %inicio contadores
    cs=0; cc=0; cf=0;
    
    % Datos iniciales
    data=[t;x];
    K=t(end);
    
    datarec=zeros(2,K);
    datarec(1,t)=t;
    datarec(2,t)=x;
    
    total_gaps=length(find(diff(t)>1));
    
    % LOCALIZADOR DE HUECOS - MODULO DE DECISIÓN de SPL para GAPS<50
    % (1 ... gap1I), (gap1F ... gap2I), (gap2F ... gap3I),...,(gapNF ... end)
    indm50=find(diff(t) > 50);
    if length(indm50)>0
        indis=[1;indm50(:)+1].';
        indfs=[indm50(:);length(t)].';
        for it=1:length(indis)
            tint=[data(1,indis(it)):data(1,indfs(it))];
            tdat = data(1,indis(it):indfs(it));
            xdat = data(2,indis(it):indfs(it));
            % RECONSTRUCCION POR SPL
            xint=interp1(tdat,xdat,tint,'*spline');
            datarec(1,tint)=tint;
            datarec(2,tint)=xint;
            cs=cs+1;
        end
        % elimino el offset inicial
        if t(1)>1, datarec(:,1:t(1)-1)=[]; end
    else
        tint=[1:K];
        xint=interp1(data(1,:),data(2,:),tint,'*spline');
        datarec(1,tint)=tint;
        datarec(2,tint)=xint;
    end
    
    
    % LOCALIZADOR DE HUECOS de CUB o FOU
    % (1 ... gap1I), (gap1F ... gap2I), (gap2F ... gap3I),...,(gapNF ... end)
    zeroo=find(datarec(1,:)==0);
    datarec(:,zeroo)=[];
    datarec2=zeros(2,K);
    datarec2(1,datarec(1,:))=datarec(1,:);
    datarec2(2,datarec(1,:))=datarec(2,:);
    indM50=find(diff(datarec(1,:)) > 50);
    if length(indM50)>0
        indIs=[1;indM50(:)+1].';
        indFs=[indM50(:);length(datarec)].';
        
        for it2=1:length(indIs)-1
            rango=300;
            gapA=indIs(it2);
            gapI=indFs(it2);
            gapF=indIs(it2+1);
            gapZ=indFs(it2+1);
            if gapI-gapA+1 > rango, gapA=gapI-rango; end
            if gapZ-gapF+1 > rango, gapZ=gapF+rango; end
            Igap.t=datarec(1, gapA:gapI);
            Igap.x=datarec(2, gapA:gapI);
            Dgap.t=datarec(1, gapF:gapZ);
            Dgap.x=datarec(2, gapF:gapZ);
            T2=[datarec(1,gapA):datarec(1,gapZ)];
            
            % CARACTERIZACIÓN, DECISION y RECONSTRUCCION por CUB o FOU
            [X2,m]=repara(Igap,Dgap,T2,p);

            if strcmp(m,'FOU')
                cf=cf+1;
            else
                cc=cc+1;
            end
            % insertar la parte reconstruida
            datarec2(1,T2)=T2;
            datarec2(2,T2)=X2;
        end
        % eliminar el offset inicial
        if datarec(1,1)>1, datarec2(:,1:datarec(1,1)-1)=[]; end
    end
    secs=toc;
    
    % guardar datos resumen
    resumen.cs=total_gaps-cc-cf;
    resumen.cc=cc;
    resumen.cf=cf;
    resumen.s=secs;
    
    % salida de la función
    tr=datarec2(1,:);
    xr=datarec2(2,:);
