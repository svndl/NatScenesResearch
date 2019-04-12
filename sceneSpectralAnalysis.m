function sceneSpectralAnalysis
    %% define path for original scenes/scrambles
    close all;
    path_dataRoot =  fullfile('/Users', 'babylab', 'Documents', 'MATLAB', '2D3D');
    path_dataScrambles = fullfile(path_dataRoot, 'scrambled');
    path_dataScenes = fullfile(path_dataRoot, 'scenes');
    %% list dir
    
    list_scenes = list_folder(path_dataScenes);
    list_scrambles = list_folder(path_dataScrambles);   
    nScenes = numel(list_scenes);
    nOri = 361;
    nFreq = 101;
    oriDiff = zeros(nOri, nScenes);
    spatDiff = zeros(nFreq, nScenes);
    std_scene = zeros(nScenes, 1);
    std_scrambled = std_scene;
    for n = 1:nScenes
       
        % load data
       load(fullfile(path_dataScrambles, list_scrambles(n).name));
       load(fullfile(path_dataScenes, list_scenes(n).name));       
       
       % convert to grayscale         
       scene_gs = rgb2gray(scene.right);
       scramble_gs = rgb2gray(scrambled);
       
       std_scene(n) = calcPixelStat_std(scene_gs);
       std_scrambled(n) = calcPixelStat_std(scramble_gs);
    
       % plot asnd save
       save_path = fullfile(path_dataRoot, ['diff_' scene.name]);
       [oriDiff(:, n), spatDiff(:, n)] = plotSpectrumDifference(scene_gs, scramble_gs, save_path, scene.name);
    end
    
    avgOriDiff = mean(oriDiff, 2);
    avgSpatDiff = mean(spatDiff, 2);
    
    figure;
    of = subplot(2, 1, 1);
    sf = subplot(2, 1, 2);

    oris = linspace(-180, 180, 361);
    frqs = linspace(0,size(scene_gs, 1)/2, 101);
        
    %% Orientation energy
    plot_options = {'Marker', 'o', 'LineStyle','none', 'Markersize', 7};
    
    plot(of, oris, avgOriDiff, plot_options{:}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'y');     
    plot(sf, frqs, avgSpatDiff, plot_options{:}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'm');  
    
    xlabel(of, 'Orientation (deg)','fontsize', 20);
    xlabel(sf, 'Cycles per Image','fontsize', 20);
    
    legend(of, {'Average Orientation Difference, %'});
    legend(sf, {'Average Spatial Difference, %'});

    set(of, 'fontsize', 20, 'box', 'off', 'LineWidth', 2);
    set(sf, 'fontsize', 20, 'box', 'off', 'LineWidth', 2);    

    title(of, 'Average Orientation Energy Difference', 'fontsize', 30, 'Interpreter', 'none');
    title(sf, 'Average Spatial Energy Difference', 'fontsize', 30, 'Interpreter', 'none');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);    
    saveas(gcf, fullfile(path_dataRoot, 'diff_Average'), 'epsc');
    close(gcf);
    
    %% plot standard dev
    figure;
    plot_options = {'Marker', 'o', 'LineStyle','none', 'Markersize', 7};
    plot(std_scene, std_scrambled, plot_options{:}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'm');
    xlabel('Scene pixel standard dev','fontsize', 20);
    ylabel('Scrambled pixel standard dev','fontsize', 20);
    set(gca, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);    
    set(gca, 'fontsize', 20, 'box', 'off', 'LineWidth', 2);
    axis square;
    title('Standard Deviation Distribution', 'fontsize', 30, 'Interpreter', 'none');
    saveas(gcf, fullfile(path_dataRoot, 'pixel_std'), 'epsc');
    close(gcf);
end
function [diffOrientation_pct, diffSpatial_pct] = plotSpectrumDifference(im1, im2, save_path, scene_name)
    
    figure('visible','on');
    oriFig = subplot(2, 2, 1);
    diffOriFig= subplot(2, 2, 2);
    sfFig = subplot(2, 2, 3);
    diffSfFig = subplot(2, 2, 4);
    
    std_im1 = calcPixelStat_std(im1);
    std_im2 = calcPixelStat_std(im2);
    
    stdLeg_1 = sprintf('Original pixel std %.3f', std_im1);
    stdLeg_2 = sprintf('Scrambled pixel std %.3f', std_im2);
    
    rawOriPow_im1 = meanOrientPow(fft2(double(im1)));
    rawOriPow_im2 = meanOrientPow(fft2(double(im2)));
    rawOriPow_diff = abs(rawOriPow_im1 - rawOriPow_im2);
    
    rawSpatPow_im1 = meanSpatialPow(fft2(double(im1)));    
    rawSpatPow_im2 = meanSpatialPow(fft2(double(im2)));
    rawSpatPow_diff = abs(rawSpatPow_im1 - rawSpatPow_im2);
    
    
    yScale = linspace(-180, 180, 361);
        
    %% Orientation energy
    plot_options = {'Marker', 'o', 'LineStyle','none', 'Markersize', 7};
    semilogy(oriFig, yScale, rawOriPow_im1, plot_options{:}, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');  
    hold(oriFig, 'on');      
    semilogy(oriFig, yScale, rawOriPow_im2, plot_options{:}, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');  
    hold(oriFig, 'on');      
    semilogy(oriFig, yScale, rawOriPow_diff, plot_options{:}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); 
    hold(oriFig, 'on');    

    xlabel(oriFig, 'Orientation (deg)','fontsize', 20);
    legend(oriFig, {stdLeg_1, stdLeg_2, 'Difference'});
    set(oriFig, 'fontsize', 20, 'box', 'off', 'LineWidth', 2);    
    title(oriFig, [scene_name ': Orientation Energy'],'fontsize', 30, 'Interpreter', 'none');
    
    %% Orientation difference plot
    diffOrientation_pct = 100*(rawOriPow_diff./rawOriPow_im1);    
    
    axis(diffOriFig);
    plot(diffOriFig, yScale, diffOrientation_pct, '-k', 'linewidth', 2);
    legend(diffOriFig, 'Scaled Orientation Difference');
    set(diffOriFig, 'Fontsize', 20, 'box', 'off', 'LineWidth', 2);
    xlabel(diffOriFig, 'Orientation (deg)','fontsize', 20);
    title(diffOriFig, 'Difference Orientation Energy, %', 'fontsize', 30);        
    
    %% Spatial energy
    frqs = linspace(0,size(im1, 1)/2, 101);
    semilogy(sfFig, frqs, rawSpatPow_im1, plot_options{:}, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');    
    hold(sfFig, 'on');    

    semilogy(sfFig, frqs, rawSpatPow_im2, plot_options{:}, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');     
    hold(sfFig, 'on');    
    
    semilogy(sfFig, frqs, rawSpatPow_diff, plot_options{:}, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    hold(sfFig, 'on');    
    
    xlabel(sfFig, 'Cycles per Image','fontsize', 20);
    set(sfFig, 'fontsize', 20, 'box', 'off', 'xlim', [-5,365],'LineWidth', 2);
    legend(sfFig, {stdLeg_1, stdLeg_2, 'Difference'});
    title(sfFig, [scene_name ': Spatial Frequency Power'], 'fontsize', 30, 'Interpreter', 'none');
    
    %% Spatial difference plot
    axis(diffSfFig);
    diffSpatial_pct = 100*(rawSpatPow_diff./rawSpatPow_im1);
    plot(diffSfFig, frqs, diffSpatial_pct, '-k', 'linewidth', 2);
    legend(diffSfFig, 'Scaled Spatial Difference');
    xlabel(diffSfFig, 'Cycles per Image','fontsize', 20);
    set(diffSfFig, 'Fontsize', 20, 'box', 'off', 'xlim', [-5,365], 'LineWidth', 2);
    title('Difference Spatial Frequency Power, %', 'fontsize', 30);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);    
    saveas(gcf, save_path, 'epsc');
    close(gcf);
end
function std_img = calcPixelStat_std(im)
    % total mean
    nPixels = numel(im);
    mean_img = sum(im(:))/(nPixels);
    
    delta_img = (im - mean_img).^2;
    sum_img = sum(delta_img(:));
    var_img = sum_img/(nPixels - 1);
    std_img = sqrt(var_img);
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
