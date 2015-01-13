%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Alejandro Santa-Cruz                                            %
%   Reconstrucción de señales de MoCap por análisis de Fourier.     %
%   Instituto Tecnológico de Buenos Aires 2011-2012                 %
%   fichero entrenamiento.m                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entrenamiento del algoritmo para obtener parámetros               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('cod');    % código
if not(isdir('./TEST')),mkdir('TEST');end

% OBTENER DATOS
addpath('data');    % contiene ficheros pts con los marcadores
fd='xbreak'; %
load xbreak.pts
z = xbreak';
% fd = 'xmaxi10';
% load xmaxi10.pts
% z = xmaxi10';
NUM_MARKERS = 34;

% [#signals   #frames]
[I K] = size(z);
% extraer los tiempos
to = z(1,:);

% Criterios de selección de marcadores
% uno aleatorio
rndm = randi(34,1,5);
% separado por bloques
trnk=[1 4 6 10];
upper=[11 14 17 22];
down=[23 26 27 30 32 34];
% separado por extremidades
hip=[1:4];
chest=[5:7];
head=[8:10];
arms=[11 12 13 17 18 19];
hands=[14 15 16 20 21 22];
legs=[23 24 25 26 29 30 31 32];
feet=[27 28 33 34];

% Selección de marcadores
work=horzcat(trnk,upper,down);
lw=length(work);

% GENERAR PRUEBA
% Numero de pruebas de la batería: P=Its*lw
Its = 100;
% modo verbose, sacar graficos, etc..
v=0;
% Matriz para almacenar los resultados
results = zeros(Its*lw,18);
% Preparacion de los ficheros de salida, escribir cabeceras
fid  = fopen([fd,'_Test_',num2str(Its*lw),'.csv'],'w');
fprintf(fid, 'gapL\tIhay\tDhay\tIamp\tDamp\tIvar\tDvar\tImad\tDmad\tMAE_F\tMAE_C\tObs\n');
fid1 = fopen([fd,'_MAEs_',num2str(Its*lw),'.csv'],'w');
fprintf(fid1, 'M\tFMae\tCMae\tSMae\tFRmse\tCRmse\tSRmse\n');

j=0;    % indice j para matriz resultado
for mrkr=1:NUM_MARKERS,
    if find(mrkr==work)     % Para cada marcador que cumpla el criterio
        fprintf(['\nmarker=',int2str(mrkr),' Observacion:\t\tCaracterizacion:\n'])
        % Extrae la serie de datos de ese marcador
        xo = z(1+mrkr,:);         % ficheros pts en una dimension [35 6000] = size(z);
        % xo = z(2+(mrkr-1)*3,:)';  % ficheros pts en tres dimensiones [103 6000] = size(z);
        
        for n=1:Its,        % Y realiza Its repeticiones
            i=j*Its + n;    % indice i para matriz resultado
            
            % ALGORITMO PROCESAMIENTO
            [data,mae,rmse,obs,gapL]=procesa(to,xo,v);
            
            % Escribir el fichero csv que se usa para calcular los parametros desde Calc Open Office
            fprintf(fid, '%d\t%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%s\n',...
                gapL,data.iay,data.day,data.idst,data.ddst,data.iv,data.dv,data.im,data.dm,(mae.fou/gapL),(mae.cub/gapL),obs.str);   
                
            % Recepcion de d de salida y almacenamientos
            results(i,1)= mrkr;     % marcador
            results(i,2)= n;        % iteracion: (mrkr,n) crean dupla IDprueba
            results(i,3)= gapL;     % tamaño gap     
            results(i,4)=data.iay;   % hay freq?
            results(i,5)=data.day;   %   "
            results(i,6)=data.idst;  % dist max
            results(i,7)=data.ddst;  %   "
            results(i,8)=data.iv;    % varianza
            results(i,9)=data.dv;    %   "
            results(i,10)=data.im;   % desviacion de la mediana
            results(i,11)=data.dm;   %   "
            results(i,12)=obs.num;  % valor observado: 1 FOU 2 CUB 3 SPL
            results(i,13)=mae.fou;  % suma de absolutos
            results(i,14)=mae.cub;  %   "
            results(i,15)=mae.spl;  %   "            
            results(i,16)=rmse.fou;  % suma de cuadrados
            results(i,17)=rmse.cub;  %   "
            results(i,18)=rmse.spl;  %   "
        end
        j=j+1;
        
        % Resultados de los MAE y RMSE para el conjunto total de P repeticiones
        tot_gaps=sum(results(:,3));
        Fmae=sum(results(:,13))./tot_gaps; Frmse=sqrt(sum(results(:,16))./tot_gaps);
        Cmae=sum(results(:,14))./tot_gaps; Crmse=sqrt(sum(results(:,17))./tot_gaps);
        Smae=sum(results(:,15))./tot_gaps; Srmse=sqrt(sum(results(:,18))./tot_gaps);
        
         % Escribir el fichero csv
        fprintf(fid1, '%d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\n',mrkr,Fmae,Cmae,Smae,Frmse,Crmse,Srmse);
        
    else
        disp(['marker=',int2str(mrkr),' no se procesa...'])
    end
end

% Cerrar y guardar prueba
fclose(fid);
fclose(fid1);
% save results

