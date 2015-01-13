%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alejandro Santa-Cruz                                                %
%   Reconstrucción de señales de MoCap por análisis de Fourier.         %
%   Instituto Tecnológico de Buenos Aires 2011-2012                     %
%   fichero ejemplo_uso_reconstruccion.m                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EJEMPLO DE USO 
%   Llamada al algoritmo que reconstruye los huecos en una señal MoCap  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('cod');
if not(isdir('IMG')),mkdir('IMG');end

% SELECCIONAR MARCADOR
addpath('data');
% load gmaxi1.pts
% z = gmaxi1';
% load gbreak16.pts
% z = gbreak16';
load gmaxi3.pts
z = gmaxi3';

% Parametros encontrados en el entrenamiento
% p=[0,0,0,0];
% cronometro
ttemp=0;

% seleccionar serie de datos del marcador
% mrkr = randi(34,1,1);
% O reconstruir todos los marcadores
for mrkr=1:34

    t = z(1,:);
    x = z(1+mrkr,:);    % si [35 6000] = size(z);        
    %x = z(2+(mrkr-1)*3,:); % si [103 6000] = size(z);

    % LLAMADA A LA RECONSTRUCCIÓN DE LA SEÑAL
    % Donde la serie x(t) es una serie de datos con huecos. Opcional params
    [tr,xr,res]=mocap_reconstr(t,x);%,p);
    
    
    % MOSTRAR RESULTADOS
    %
    disp(['Marker[', num2str(mrkr),'] Metodos  SPL:', int2str(res.cs), ' CUB:',int2str(res.cc), ' FOU:',int2str(res.cf),...
    '     t=',num2str(res.s),' secs'])
    ttemp = ttemp + res.s;
    name = ['./IMG/rec_marker',num2str(mrkr)];
    
    % Guardar fichero
    dlmwrite([name,'.txt'] ,xr, '\n');
    % Guardar gráfico
    old=zeros(1,t(end));
    old(1,t)=t;
    old(2,t)=x;
    ind=find(old(1,:)==0);
    old(:,ind)=NaN;
    figure(1)
    set(gcf,'Visible','off');
    clf
    subplot(211)
    title('Secuencia reconstruida')
    plot(tr,xr,'g')
    hold on
    plot(old(1,:),old(2,:),'k')
    xlabel('Tiempo (frames)');
    ylabel('Posicion');
    legend('Reconstr','Original','Location','NorthEast');
    print(gcf,'-dpng','-r125',[name,'.png']);
end

disp(['Tiempo total reconstrucciones:   ',num2str(ttemp), ' secs'])
