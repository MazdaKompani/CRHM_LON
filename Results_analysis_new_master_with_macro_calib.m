
% INPUT DATA -----------------------------------

%In lines 32 and 62 we can change the column number and the name of the
%simulated varaible
% in Line     we can change the onbserved value
% load observations

climatedata = importdata('LON2011-2018-V3-1-climate.obs');
time = datenum(climatedata.data(:,1) + 693960);
precip = timeseries(climatedata.data(:,6),time);
temp = timeseries(climatedata.data(:,2),time);

Qsurfdata = importdata('LON2011-2018-V3-1-Qsurf.obs');
time_QSurf = datenum(Qsurfdata.data(:,1) + 693960);
QSurf = timeseries(Qsurfdata.data(:,2),time_QSurf);

QTiledata = importdata('LON2011-2018-V3-1-Qtile.obs');
time_QTile = datenum(QTiledata.data(:,1) + 693960);
QTile = timeseries(QTiledata.data(:,2),time_QTile);

GWLdata = importdata('LON2011-2018-V3-1-GWL.obs');
time_GWL = datenum(GWLdata.data(:,1) + 693960);
GWL = timeseries(GWLdata.data(:,2),time_GWL);

porosity = 0.0333;

% load CRHM results

%model_tests_all = {'Lattest version','CRHM_output_1.txt';...
%                'initial (Mazda) + parametersUpdated','CRHM_output_2.txt';...;
%                'initial (Mazda) + removed macro','CRHM_output_3.txt';...
%                'initial (Mazda) + removed macro + parameters updated','CRHM_output_4.txt';...
%                'Old model -> new tillage module','CRHM_output.txt'};
            
model_tests_all = {'Lattest version','CRHM_output_1.txt'};
            
model_tests_2plot = [1];

model_tests = model_tests_all(model_tests_2plot,:);
            
crhm_col_QSurf = 3; % col including time
crhm_col_QSTile = 17;% % col including time
crhm_col_GWL = 4; %% col including time
crhm_col_evap = 9;

% -------------------------------------------------------------

% Plot results
for i=1:numel(model_tests(:,1))
   figure('name',['Model ',num2str(i)])
   h1 = axes;
    plot(QSurf,'ro','Markersize',2)
    hold on
    
    CRHMdata_model = importdata(model_tests{i,2}); % initial (Mazda)
    
    time_chrm = datenum(CRHMdata_model.data(:,1) + 693960);
    QSurf_crhm_i = timeseries(CRHMdata_model.data(:,crhm_col_QSurf),time_chrm);
    QTile_crhm_i = timeseries(CRHMdata_model.data(:,crhm_col_QSTile),time_chrm);
    GWL_crhm_i = timeseries(CRHMdata_model.data(:,crhm_col_GWL),time_chrm);
    crhm_hru_actet_i = timeseries(CRHMdata_model.data(:,crhm_col_evap),time_chrm);
  
    plot(QSurf_crhm_i,'k')
    ax1 = gcf;
    hold on
   
    ylabel('Surface flow (mm/h)')
    %ylabel('[mm/h]')
    ylim([0 30])
    xlim([time_chrm(1) time_chrm(end)])
    grid on
    datetick('x','mmm-yyyy','keeplimits','keepticks')
    
    legend('Obs. ',...
    ['Simul. ']) 
    
    % plot precip
    h2 = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2, 'Ydir', 'reverse')
    set(h2, 'YAxisLocation', 'Right')
    set(h2, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2, 'XLim', get(h1, 'XLim'))
    set(h2, 'YLim', [-30 150])
    set(h2, 'Color', 'None')
    set(h2, 'Xtick', [])
    alpha 0.1

    legend('Temperature (^oC)',...
        'Precipitation (mm)')
    
    
    
    
    %title ({['EOF runoff - Model ',num2str(i)],model_tests{i,1}})
    
    % Performence calculation
     [Nash,RMSE,Bias,dataUse,CHRMout_data_interp] = PerformenceCalculator(QSurf_crhm_i.Time,QSurf_crhm_i.Data,...
         QSurf.Time,QSurf.Data);
     PerfText = {'Performence:',...
                ['NSE =', num2str(Nash)],...
                ['RMSE =', num2str(RMSE)],...
                ['Bias =', num2str(Bias)]};
     dim = [.2 .5 .3 .3];
     annotation('textbox',dim,'String',PerfText,'FitBoxToText','on','backgroundColor','w');
    
     saveas(gcf,['Results_runoff_',model_tests{i,1}],'png')
     
     
    %%%% figure;
     
    %%%% scatter(dataUse,CHRMout_data_interp)
    %%%% title('Surface flow (mm/h)')
     
end

%plotting time series of observed and simulated surface flow _______

figure
hs = axes;
plot(QSurf.Time,QSurf.Data,'ro','Markersize',2)
hold on
% %plot(QSurf_crhm_1,'k')
% hold on
% %plot(QSurf_crhm_2,'g:')
% hold on
plot(QSurf_crhm_i,'k','Markersize',2)

grid on
ylabel('Surface flow [mm/h]')
%title('Surface flow (mm/h)')
set (hs,'XLim',[datenum([2014,11,5,1,1,0]), datenum([2016,4,5,1,1,1])])
%xlim([datenum([2015,4,5,1,1,0]), datenum([2016,4,5,1,1,1])])
datetick('x','mmm-yyyy', 'keeplimits')
%% Performence calculation
%     [NashS,RMSES,BiasS,dataUseS,CHRMout_data_interpS] = PerformenceCalculator(QSurf_crhm_i.Time,QSurf_crhm_i.Data,...
%         QSurf.Time,QSurf.Data);
%     PerfTextS = {'Performence:',...
%                ['NSE =', num2str(NashS)],...
%                ['RMSE =', num2str(RMSES)],...
%                ['Bias =', num2str(BiasS)]};
%     dimS = [.2 .5 .3 .3];
%     annotation('textbox',dimS,'String',PerfTextS,'FitBoxToText','on','backgroundColor','w');
    
legend('Obs. ',...
    ['Simul. ']) 

     % plot precip
    h2 = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2, 'Ydir', 'reverse')
    set(h2, 'YAxisLocation', 'Right')
    set(h2, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2, 'XLim', get(hs, 'XLim'))
    set(h2, 'YLim', [-30 150])
    set(h2, 'Color', 'None')
    set(h2, 'Xtick', [])
    alpha 0.1

    legend('Temperature (^oC)',...
        'Precipitation (mm)')
    




%plotting time series of observed and simulated TILE flow________________
%_
%legend('Obs','Model_1','Model_2','Model_3')

figure
hT=axes;
plot(QTile.Time,QTile.Data,'ro','Markersize',2)
hold on
% %plot(QSurf_crhm_1,'k')
% hold on
% %plot(QSurf_crhm_2,'g:')
% hold on
plot(QTile_crhm_i,'k')
datetick('x','mmm-yyyy')
grid on
ylabel('Tile flow [mm/h]')
%title('Tile flow (mm/h)')

set(hT, 'YLim', [0 0.7])
set(hT, 'XLim', get(h1, 'XLim'))
%xlim([time_chrm(1) time_chrm(end)])    
% Performence calculation
     [NashT,RMSET,BiasT,dataUseT,CHRMout_data_interpT] = PerformenceCalculator(QTile_crhm_i.Time,QTile_crhm_i.Data,...
         QTile.Time,QTile.Data);
     PerfTextT = {'Performence:',...
                ['NSE =', num2str(NashT)],...
                ['RMSE =', num2str(RMSET)],...
                ['Bias =', num2str(BiasT)]};
     dimT = [.2 .5 .3 .3];
     annotation('textbox',dimT,'String',PerfTextT,'FitBoxToText','on','backgroundColor','w');

     legend('Obs. ',...
    ['Simul. ']) 
     % plot precip
    h2 = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2, 'Ydir', 'reverse')
    set(h2, 'YAxisLocation', 'Right')
    set(h2, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2, 'XLim', get(h1, 'XLim'))
    set(h2, 'YLim', [-30 150])
    set(h2, 'Color', 'None')
    set(h2, 'Xtick', [])
    alpha 0.1

   legend('Temperature (^oC)',...
        'Precipitation (mm)')

     
     

% legend('Obs',... %'Model_1','Model_2',.. 
% 'Model_3')
% 
% %plot(QTile,'r.')
% %datetick('x','mmm-yyyy')
% %grid on
% %ylabel('[mm/h]')
% %title('Tile flow (mm/h)')
% %legend('Obs','Model_1','Model_2','Model_3')
% 

%________________________TILE flow plot partly

figure
hTT=axes;
plot(QTile.Time,QTile.Data,'ro','Markersize',2)
hold on
% %plot(QSurf_crhm_1,'k')
% hold on
% %plot(QSurf_crhm_2,'g:')
% hold on
plot(QTile_crhm_i,'k')
datetick('x','mmm-yyyy')
grid on
ylabel('Tile flow [mm/h]')
%title('Tile flow (mm/h)')

set(hTT, 'YLim', [0 0.7])
set(hTT, 'XLim', [datenum([2014,11,5,1,1,0]), datenum([2016,4,5,1,1,1])])
%xlim([time_chrm(1) time_chrm(end)])    
%datetick('x','mmm-yyyy', 'keeplimits','keepticks')
datetick('x','mmm-yyyy','keeplimits')
     
legend('Obs. ',...
    ['Simul. ']) 
     % plot precip
    h2TT = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2TT, 'Ydir', 'reverse')
    set(h2TT, 'YAxisLocation', 'Right')
    set(h2TT, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2TT, 'XLim', get(hTT, 'XLim'))
    set(h2TT, 'YLim', [-30 150])
    set(h2TT, 'Color', 'None')
    set(h2TT, 'Xtick', [])
    alpha 0.1

   legend('Temperature (^oC)',...
        'Precipitation (mm)')

     
     







%_____________________plotting GWL
figure
%subplot(311)
hG=axes;
plot(GWL,'r.')
hold on
GWL_crhm_i_conv_to_m = GWL_crhm_i / porosity / 1000 - (3-1.5); % well checked and correct (Mazda and Diogo)
plot(GWL_crhm_i_conv_to_m,'k')
%hold on
% plot(GWL_crhm_2,'g:')
% hold on
% plot(GWL_crhm_3,'r:')
datetick('x','mmm-yyyy')
grid on
ylabel('Groundwater/soil water level [m]')
%title('GW level (m)')
set(hG, 'YLim', [-1.5 4.7])
set(hG, 'XLim', get(h1, 'XLim'))

% Performence calculation
     [NashG,RMSEG,BiasG,dataUseG,CHRMout_data_interpG] = PerformenceCalculator(GWL_crhm_i.Time,GWL_crhm_i_conv_to_m.Data,...
         GWL.Time,GWL.Data);
     PerfTextG = {'Performence:',...
                ['NSE =', num2str(NashG)],...
                ['RMSE =', num2str(RMSEG)],...
                ['Bias =', num2str(BiasG)]};
     dimG = [.2 .5 .3 .3];
     annotation('textbox',dimG,'String',PerfTextG,'FitBoxToText','on','backgroundColor','w');

% legend('Obs','Model_1','Model_2','Model_3')
% 
% plot(GWL,'r.')
% datetick('x','mmm-yyyy')
% grid on
% ylabel('[m]')
% title('GW level (m)')
legend('Obs. ',...
    ['Simul. ']) 
  
     % plot precip
    h2 = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2, 'Ydir', 'reverse')
    set(h2, 'YAxisLocation', 'Right')
    set(h2, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2, 'XLim', get(h1, 'XLim'))
    set(h2, 'YLim', [-30 150])
    set(h2, 'Color', 'None')
    set(h2, 'Xtick', [])
    alpha 0.1

   legend('Temperature (^oC)',...
        'Precipitation (mm)')

%_____________________plotting GWL    partly_____________________
figure
%subplot(311)
hGG=axes;
plot(GWL,'r.')
hold on
GWL_crhm_i_conv_to_m = GWL_crhm_i / porosity / 1000 - (3-1.5); % well checked and correct (Mazda and Diogo)
plot(GWL_crhm_i_conv_to_m,'k')
%hold on
% plot(GWL_crhm_2,'g:')
% hold on
% plot(GWL_crhm_3,'r:')
datetick('x','mmm-yyyy')
grid on
ylabel('Groundwater/soil water level [m]')
%title('GW level (m)')
set(hGG, 'YLim', [-1.5 4.7])
set(hGG, 'XLim', [datenum([2014,11,5,1,1,0]), datenum([2016,4,5,1,1,1])])

datetick('x','mmm-yyyy','keeplimits')
legend('Obs. ',...
    ['Simul. ']) 
  
     % plot precip
    h2GG = axes;
    plot(temp,'Color',[0.7 0.7 0.7],'marker','.')
    hold on
    plot(precip,'Color','b','linewidth',2)
    set(h2GG, 'Ydir', 'reverse')
    set(h2GG, 'YAxisLocation', 'Right')
    set(h2GG, 'XAxisLocation', 'top')
    ylabel('Precipitation (mm) and Temperture (^oC)')
    
    set(h2GG, 'XLim', get(hGG, 'XLim'))
    set(h2GG, 'YLim', [-30 150])
    set(h2GG, 'Color', 'None')
    set(h2GG, 'Xtick', [])
    alpha 0.1

   legend('Temperature (^oC)',...
        'Precipitation (mm)')

     
     
%__________________________plotting the cumulative curves

figure;
subplot(131)
plot(cumsum(precip.Data))
hold on
plot(cumsum(QTile_crhm_i.Data))
hold on
plot(cumsum(crhm_hru_actet_i.Data))
hold on; 
plot(cumsum(QSurf_crhm_i.Data))
legend('precip','QTile_crhm_i','crhm_hru_actet_i','QSurf_crhm_i')
grid on
subplot(132)
CQTile_crhm_i.Data=cumsum(QTile_crhm_i.Data)      %added by Mazda
plot(cumsum(QTile_crhm_i.Data))
hold on
CQTile.Data=cumsum(QTile.Data)      %added by Mazda
plot(cumsum(QTile.Data))
title('Tile flow')
legend('Simul.','obs.')
grid on
subplot(133)
CQSurf_crhm_i.Data=cumsum(QSurf_crhm_i.Data)      %added by Mazda
plot(cumsum(QSurf_crhm_i.Data))
hold on
CQSurf.Data=cumsum(QSurf.Data)                   %added by Mazda
plot(cumsum(QSurf.Data))
grid on
title('Surface runoff')
legend('Simul.','obs.')


% Total Performence calculation
     [NashTo,RMSETo,BiasTo,dataUseTo,CHRMout_data_interpTo] = PerformenceCalculatorTo(QSurf_crhm_i.Time,QSurf_crhm_i.Data,...
         QSurf.Time,QSurf.Data,QTile_crhm_i.Time,QTile_crhm_i.Data,QTile.Time,QTile.Data,GWL_crhm_i.Time,...
         GWL_crhm_i_conv_to_m.Data,GWL.Time,GWL.Data);
     PerfTextTo = {'Performence:',...
                ['NSE =', num2str(NashTo)],...
                ['RMSE =', num2str(RMSETo)],...
                ['Bias =', num2str(BiasTo)]};
     dimTo = [.2 .5 .3 .3];
     annotation('textbox',dimTo,'String',PerfTextTo,'FitBoxToText','on','backgroundColor','w');

     
     % Performence calculation
     [NashCT,RMSECT,BiasCT,dataUseCT,CHRMout_data_interpCT] = PerformenceCalculator(QTile_crhm_i.Time,CQTile_crhm_i.Data,...
         QTile.Time,CQTile.Data);
     PerfTextT = {'Cumul. Tile flow performance:',...
                ['NSE =', num2str(NashCT)],...
                ['RMSE =', num2str(RMSECT)],...
                ['Bias =', num2str(BiasCT)]};
     dimT = [.2 .5 .3 .3];
     annotation('textbox',dimT,'String',PerfTextT,'FitBoxToText','on','backgroundColor','w');

     
     
     % Performence calculation
     [NashCS,RMSECS,BiasCS,dataUseCS,CHRMout_data_interpCS] = PerformenceCalculator(QSurf_crhm_i.Time,CQSurf_crhm_i.Data,...
         QSurf.Time,CQSurf.Data);
     PerfTextT = {'Cumul. Surface flow Performence:',...
                ['NSE =', num2str(NashCS)],...
                ['RMSE =', num2str(RMSECS)],...
                ['Bias =', num2str(BiasCS)]};
     dimT = [.2 .5 .3 .3];
     annotation('textbox',dimT,'String',PerfTextT,'FitBoxToText','on','backgroundColor','w');

     

% Peformence Calculator (Nash, RMSE, Bias)
function [Nash,RMSE,Bias,dataUse,CHRMout_data_interp] = PerformenceCalculator(CRHM_time,CRHMout_data_var,time,data)

CHRMout_data_interp = interp1(datenum(CRHM_time),CRHMout_data_var,datenum(time));
% remove NaNs
isnanLoc = isnan(CHRMout_data_interp);
CHRMout_data_interp(isnanLoc) = [];
dataUse = data;
dataUse(isnanLoc) = [];

% Calculate Performence
% Nash
numerator=(dataUse-CHRMout_data_interp).^2;
denominator=(data-mean(data)).^2;
Nash =1-(sum(numerator)/sum(denominator));
% RMSE
Sumcal = (CHRMout_data_interp-dataUse).^2;
numerator = sum(Sumcal);
n=numel(dataUse);
RMSE=(numerator/n)^(1/2);
% Bias
numerator = sum(dataUse);
denominator = sum(CHRMout_data_interp);
Bias = numerator/denominator-1;

end

% Peformence Calculator (Nash, RMSE, Bias) for tile, surface flows and GWL
function [NashTo,RMSETo,BiasTo,dataUseTo_No,CHRMout_data_interpTo_No] = PerformenceCalculatorTo(CRHM_timeS,CRHMout_data_varS,...
    timeS,dataS,CRHM_timeT,CRHMout_data_varT,timeT,dataT,CRHM_timeG,CRHMout_data_varG,timeG,dataG)

CHRMout_data_interpS = interp1(datenum(CRHM_timeS),CRHMout_data_varS,datenum(timeS));
% remove NaNs
isnanLocS = isnan(CHRMout_data_interpS);
CHRMout_data_interpS(isnanLocS) = [];
dataUseS = dataS;
dataUseS(isnanLocS) = [];

CHRMout_data_interpS_No=normalize(CHRMout_data_interpS)
dataUseS_No=normalize(dataUseS)
dataS_No=normalize(dataS)
NumS=numel(dataUseS_No)



CHRMout_data_interpT = interp1(datenum(CRHM_timeT),CRHMout_data_varT,datenum(timeT));
% remove NaNs
isnanLocT = isnan(CHRMout_data_interpT);
CHRMout_data_interpT(isnanLocT) = [];
dataUseT = dataT;
dataUseT(isnanLocT) = [];

CHRMout_data_interpT_No=normalize(CHRMout_data_interpT)
dataUseT_No=normalize(dataUseT)
dataT_No=normalize(dataT)
NumT=numel(dataUseT_No)


CHRMout_data_interpG = interp1(datenum(CRHM_timeG),CRHMout_data_varG,datenum(timeG));
% remove NaNs
isnanLocG = isnan(CHRMout_data_interpG);
CHRMout_data_interpG(isnanLocG) = [];
dataUseG = dataG;
dataUseG(isnanLocG) = [];

CHRMout_data_interpG_No=normalize(CHRMout_data_interpG)
dataUseG_No=normalize(dataUseG)
dataG_No=normalize(dataG)
NumG=numel(dataUseG_No)


CHRMout_data_interpTo_No(1:NumS)=CHRMout_data_interpS_No(1:NumS)
dataUseTo_No(1:NumS)=dataUseS_No(1:NumS)
dataTo_No(1:NumS)=dataS_No(1:NumS)

CHRMout_data_interpTo_No(NumS+1:NumS+NumT)=CHRMout_data_interpT_No(1:NumT)
dataUseTo_No(NumS+1:NumS+NumT)=dataUseT_No(1:NumT)
dataTo_No(NumS+1:NumS+NumT)=dataT_No(1:NumT)

CHRMout_data_interpTo_No(NumS+NumT+1:NumS+NumT+NumG)=CHRMout_data_interpG_No(1:NumG)
dataUseTo_No(NumS+NumT+1:NumS+NumT+NumG)=dataUseG_No(1:NumG)
dataTo_No(NumS+NumT+1:NumS+NumT+NumG)=dataG_No(1:NumG)





% Calculate Performence
% Nash
numeratorTo_No=(dataUseTo_No-CHRMout_data_interpTo_No).^2;
denominatorTo_No=(dataTo_No-mean(dataTo_No)).^2;
NashTo =1-(sum(numeratorTo_No)/sum(denominatorTo_No));
% RMSE
SumcalTo_No = (CHRMout_data_interpTo_No-dataUseTo_No).^2;
numeratorTo_No = sum(SumcalTo_No);
nTo_No=numel(dataUseTo_No);
RMSETo=(numeratorTo_No/nTo_No)^(1/2);
% Bias
numeratorTo_No = sum(dataUseTo_No);
denominatorTo_No = sum(CHRMout_data_interpTo_No);
BiasTo = numeratorTo_No/denominatorTo_No-1;

end
