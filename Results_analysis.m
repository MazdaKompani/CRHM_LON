
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

porosity = 0.15;

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
  
    plot(QSurf_crhm_i,'k')
    ax1 = gcf;
    hold on
   
    ylabel('Obs and Simulated surface runoff (mm/h)')
    %ylabel('[mm/h]')
    ylim([0 30])
    xlim([time_chrm(1) time_chrm(end)])
    grid on
    datetick('x','mmm-yyyy','keeplimits','keepticks')
    
    legend('Obs (runoff, left axis, mm/h)',...
    ['Model_',num2str(i), 'runoutflow-m, left axis],mm/h)'])
    
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
    set(h2, 'YLim', [-30 60])
    set(h2, 'Color', 'None')
    set(h2, 'Xtick', [])
    alpha 0.1

    legend('Precipitation (right axis, mm)',...
        'Temperature (right axis, ^oC)')
    
    title ({['EOF runoff - Model ',num2str(i)],model_tests{i,1}})
    
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
     
     
     figure;
     scatter(dataUse,CHRMout_data_interp)
     
end

%legend('Obs','Model_1','Model_2','Model_3')

figure
plot(QTile,'ro')
% hold on
% %plot(QSurf_crhm_1,'k')
% hold on
% %plot(QSurf_crhm_2,'g:')
% hold on
plot(QTile_crhm_i,'k')
% datetick('x','mmm-yyyy')
% grid on
% ylabel('[mm/h]')
% title('Tile flow (mm/h)')
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
figure
%subplot(311)
plot(GWL,'r.')
hold on
GWL_crhm_i_conv_to_m = GWL_crhm_i / porosity / 1000 - 1.5; % well checked and correct (Mazda and Diogo)
plot(GWL_crhm_i_conv_to_m,'k')
% hold on
% plot(GWL_crhm_2,'g:')
% hold on
% plot(GWL_crhm_3,'r:')
% datetick('x','mmm-yyyy')
% grid on
% ylabel('[m]')
% title('GW level (m)')
% legend('Obs','Model_1','Model_2','Model_3')
% 
% plot(GWL,'r.')
% datetick('x','mmm-yyyy')
% grid on
% ylabel('[m]')
% title('GW level (m)')



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
