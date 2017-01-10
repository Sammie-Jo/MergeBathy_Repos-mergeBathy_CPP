%Create Figs
% The version of this in Test_Center is better!!!!
function createBathyCppFigs


% hold on;
% axis image
% [cmin, cmax]=caxis;
% NCEX
%  caxis([-700,0]);
%  caxis([-50,50]);

%  hold on;
%  [cmin, cmax]=caxis;
%  caxis([cmin/10,cmax/10]);

% DBDBV
%  caxis([-3000,0]);
%  caxis([-50,50]);



%%
% Directory for Test Set to run
% TestSetRoot='Test_Set_01_runValidationData_Kevins/Active_Testing_Site/';

% Directory for output files from MergeBathy to plot
% outputDir = 'output_files_win/';
% % bitFlag = 'x86/';
% bitFlag = 'x64/';
% configFlag = 'Debug';
% % configFlag = 'Release';
% loc = ([[outputDir bitFlag] configFlag]);
% 
% figDir = 'figs';
% % [s, m1, m2] = mkdir(['../Active_Testing_site/' loc], figDir);
% loc = [loc '/'];
% figDir = [figDir '/'];

filename1 = 'T01C04_CPP_DBDBV_0.5_NoOverlap';
filename2 = 'T01C04_CPP_DBDBV_0.5_NoOverlapMBZ926';
% filename2 = 'T01C04_CPP_DBDBV_0.5_NoOverlap_GMT_xyz_xyd';
% filename2 = 'T01C04_CPP_DBDBV_0.5_NoOverlap_GMT_xyz_xyde';
% filename1 = 'T01C05ehv_CPP_PORTO_185x185_MBZ400x400';
% filename2 = 'T01C05ehv_CPP_PORTO_400x400_MBZ_xyde';
% filename1 = 'T01C05ehv_CPP_PORTO_185x185_MBZ100x100';
% filename2 = 'T01C05ehv_CPP_PORTO_100x100_MBZ_xyde';

fList= char([filename1 '.txt'],...
	[filename2 '.txt']);

 bathyTitles = char([filename1 'Bathymetry'],...
	[filename2 'Bathymetry']);
 
 uncertTitles = char([filename1 'Uncertainty'],...
	[filename2 'Uncertainty']);

  bathyFigs = char('NCEX50x50_Bathy',...
	'NCEX100x100_Bathy');
 
 uncertFigs = char('NCEX50x50_Uncertainty',...
     'MTrenchUncertainty_Krig');

 
 [r,c] = size(fList);
 for cnt=1:r
    cFile = fList(cnt,:);
disp(cFile)
    %Generate Bathymetry Figures===========================================
    A = geoDataArray(cFile);
    figure; surf(geoSign(A), 1000, 1000, 'regrid'); 
    %set(gca, 'Clim', [-11000 -6000]);
%   SaveFig( bathyTitles(cnt,:), [[loc figDir] bathyFigs(cnt,:)], 'Bathymetry');
    LabelFig(bathyTitles(cnt,:),'Bathymetry');
%     SaveFig( [figDir bathyFigs(cnt,:)] );
%     FormatFig(cnt,1, bathyFigs(cnt,:));
    
    %Generate Uncertainty Figures==========================================
    B = geoDataArray(cFile, 4);
    figure; surf(B, 1000, 1000, 'regrid'); 
%    set(gca, 'Clim', [0 20]);
%    figure; plot(B, 'ColorZ')
    
%     c = bboxcut(b, 13.1, 12.9, 146, 146.2)  
% %    figure; plot(c)
% %    figure; plot(c, 'colorz')
% 
%     figure; surf(c, 1000, 1000, 'regrid')
%     figure; surf(geolimit(c, 10) , 1000, 1000, 'regrid')
%  
%     SaveFig( uncertTitles(cnt,:), [[loc figDir] uncertFigs(cnt,:)], 'Uncertainty');
    LabelFig(uncertTitles(cnt,:),'Uncertainty');
%     SaveFig( [[loc figDir] uncertFigs(cnt,:)] );
%     FormatFig(cnt,2, uncertFigs(cnt,:));
    
    if(cnt==3) 
        D1=B;
    elseif(cnt==4)
        D3=B;
    elseif(cnt==9)
        D2=B;
    elseif(cnt==10)
        D4=B;
    end
 end
 
%   % Print DBDBV uncertainty differences when kriging
%   D = D2 - D1; %DBDBV_Uncertainty - DBDBV_Uncertainty_Krig
%   figure; surf(D, 1000, 1000, 'regrid'); 
% %   SaveFig( 'DBDBV Kriging Uncertainty Differences', [[loc figDir] 'DBDBV_Kriging_Uncert_Diffs'], 'Uncertainty');
%   LabelFig('DBDBV Kriging Uncertainty Differences','Uncertainty');
% %   SaveFig( [[loc figDir] 'DBDBV_Kriging_Uncert_Diffs'] );
%   FormatFig(cnt+1,2,[[loc figDir] 'DBDBV_Kriging_Uncert_Diffs']);
% 
%   
%   D = D4 - D3; %DBDBV_Uncertainty_NoOverlap - DBDBV_Uncertainty_Krig_NoOverlap
%   figure; surf(D, 1000, 1000, 'regrid'); 
%   LabelFig('DBDBV Kriging Uncertainty No Overlap Differences','Uncertainty');
% %   SaveFig( [[loc figDir] 'DBDBV_Kriging_Uncert_NoOverlap_Diffs'] );
%   FormatFig(cnt+1,2,[[loc figDir] 'DBDBV_Kriging_Uncert_NoOverlap_Diffs']);

  
end

 function LabelFig(figTitle, yLabel)
 
  % Figure Formating
    xlabel('Longitude', 'FontSize', 14)
    ylabel('Latitude', 'FontSize', 14)
    set(gca, 'FontSize', 14)
    title(figTitle, 'FontSize', 14)
    h = colorbar;
    set(h, 'FontSize', 14)
    set(get(h, 'YLabel'), 'String', strcat(yLabel,' (m)'), 'FontSize', 14)
 end
 
 function SaveFig(figFile)

    cFigFile = figFile;
    cnt2=1;
    while (exist(strcat(cFigFile,'.jpg'),'file')~=0)
        cnt2= cnt2 + 1;
        cFigFile= strcat(figFile,'_', int2str(cnt2));
    end
    print('-djpeg', '-r300', strtrim(cFigFile));
    hgsave(strtrim(cFigFile));
 end
 
 function FormatFig(cnt,fflag,figFile)
    if((cnt==1 || cnt==2) || (cnt==7 || cnt==8)) %NCEX
        hold on;
        axis image;
        if fflag==1
            caxis([-700,0]);
        else 
            caxis([-50,50]);
        end
        hold off;
%         SaveFig( strcat(figFile, '_Formatted'));
    elseif((cnt==3 || cnt==4) || (cnt==9 || cnt==10) || cnt==13)%DBDBV
        hold on;
        axis image;
        if fflag==1
            caxis([-3000,0]);
        else
            caxis([-50,50]);
        end
        hold off;
%         SaveFig( strcat(figFile, '_Formatted'));
    end
 end
