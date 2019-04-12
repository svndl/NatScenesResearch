function toPlotMeans = freqAnalyser(raw1, raw2, saveString)
    fig_an(1) = figure('visible','on');
    oriFig(1) = subplot(2, 2, 1);
    oriFig(2) = subplot(2, 2, 2);
    sfFig(1) = subplot(2, 2, 3);
    sfFig(2) = subplot(2, 2, 4);
    
    %meanOriPow = meanOrientPow(meanAmp);
    
    rawOriPow{1} = meanOrientPow(fft2(double(raw1)));
    %filterOriPow{1} = meanOrientPow(fft2(double(filter1)));
    %maskOriPow{1} = meanOrientPow(fft2(double(mask1)));
    
    rawOriPow{2} = meanOrientPow(fft2(double(raw2)));
    %filterOriPow{2} = meanOrientPow(fft2(double(filter2)));
    %maskOriPow{2} = meanOrientPow(fft2(double(mask2)));
    
    yScale = linspace(-180, 180, 361);
    figTitles = {'scene','scramble','diff'};
    
    %% disp Orientation energy
    for z=1:2
        h_c = semilogy(oriFig(z), yScale, rawOriPow{z});
        hold on;
        title(oriFig(z), [figTitles{z},': Orientation Energy'],'fontsize',20);
        xlabel(oriFig(z), 'Orientation (deg)','fontsize',7);
        set(h_c(1),'Marker','o','LineStyle','none','markersize',5);
        set(oriFig(z),...
            'fontsize',5,'box','off', ...
            'ylim',[10^-15,10^5],...
            'tickdir','out','LineWidth',.5,'TickLength',[0.0001 0.0001],'XMinorTick','off','YMinorTick','off');
        hold off;
    end
    
    %% disp spatial energy
    
    rawSpatPow{1} = meanSpatialPow(fft2(double(raw1)));
    
    rawSpatPow{2} = meanSpatialPow(fft2(double(raw2)));
    
    frqs = linspace(0,size(raw1, 1)/2, 101);
    for z=1:2
        h_c = semilogy(sfFig(z), frqs, rawSpatPow{z});
        hold on
        set(h_c(1),'Marker','o','LineStyle','none','markersize',5);
        title(sfFig(z), [figTitles{z},': Spatial Frequency Power'],'fontsize',7);
        xlabel(sfFig(z), 'Cycles per Image','fontsize',7);
        legend(sfFig(z), 'mean', 'raw', 'filtered','masked','location','SouthEast');
        set(sfFig(z),...
            'fontsize',5,'box','off', ...
            'ylim',[10^-15,10^5],'xlim',[-5,365],...
            'tickdir','out','LineWidth',.5,'TickLength',[0.0001 0.0001],'XMinorTick','off','YMinorTick','off');        
        hold off
    end
    %% make figure 2
    fig_an(2) = figure('visible','on');
    toPlot(:, 1) = (rawSpatPow{1}-rawSpatPow{2})./rawSpatPow{1};
    toPlotMeans = mean(abs(toPlot));
    figure(fig_an(2));
    hold on
    plot(frqs, toPlot);
    title([figTitles{3},': Spatial Frequency Power'],'fontsize',7);
    xlabel('Cycles per Image','fontsize',7);
    legend(['raw ',num2str(toPlotMeans(1),'%1.3f')], ...
            'location','SouthEast');
    set(gca,...
        'fontsize',5,'box','off','xlim',[-5,365], ...
        'tickdir','out','LineWidth',.5,'TickLength',[0.0001 0.0001],'XMinorTick','off','YMinorTick','off');  
    %% save figures
    for z=1:2
        lastSlash = max(strfind(saveString,'/'));
        saveName = [saveString(lastSlash+1:lastSlash+4),num2str(z),saveString(lastSlash+5:end)];
        if exist('plot2svg','file')
            plot2svg([saveString(1:lastSlash),saveName,'.svg'],fig_an(z));
            system(['cd ',saveString(1:lastSlash), ...
                    '; svg2png ',saveName,'.svg' ...
                    '; rm ',saveName,'.svg']);
                    % note: for some reason plot2svg screws up the legend position when
                    % run in parallel with parfor
        else
            saveas(fig_an(z), [saveString(1:lastSlash),saveName], 'jpeg');
        end
        close(fig_an(z));
    end
end
    
function out = meanOrientPow(inImage)
    [nx, ny] = size(inImage);
    [x, y] = ndgrid(linspace(-1, 1, nx),linspace(-1, 1, ny));
    [th, r] = cart2pol(x, y);
    
    th = round( (180/pi)*th);
    th = th(:);
    out = zeros(361, 1);
    idx = 1;
        
    f = fftshift(inImage);
    th(r < .1) = 255;
    th(r > 1) = 255; 
    for i = -180:180          
        %normFac = max(r(th==i));
        out(idx) = mean(abs(f(th == i)));
        idx = idx + 1;
    end
end
function out = meanSpatialPow(inImage)
    [nx, ny] = size(inImage);
    [x, y] = ndgrid(linspace(-1, 1, nx), linspace(-1, 1, ny));
    [th, r] = cart2pol(x,y);
    
    r = round(100*r);
    out = zeros(101, 1);
    idx = 1;
    f = fftshift(inImage);
    for i = 0:100
        out(idx) = mean(abs(f(r == i)));
        idx = idx + 1;
    end
end
