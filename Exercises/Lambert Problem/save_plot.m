function [] = save_plot(pdf_plot, filename, dpi)

    % save_plot.m - Astronomical constants of celestial bodies.
    %
    % PROTOTYPE:
    %	[] = save_plot(pdf_plot, filename, dpi)
    %
    % DESCRIPTION:
    %   This function saves a plot as a fig, pdf and png file.
    %
    % INPUT:
    %   pdf_plot[1]     Figure handle
    %   filename        Name of the file in string
    %   dpi[1]          Dots per inch (dpi) of the saved file
    %
    % OUTPUT:
    %   Figures
    %
    % CALLED FUNCTIONS:
    %   None
    %
    % AUTHOR:
    %   Yi Qiang Ji Zhang, 23/03/2023, MATLAB, save_plot.m
    %
    % -------------------------------------------------------------------------

    % Save figure
    savefig([filename, '.fig'])

    % Save pdf
    set(pdf_plot, 'Units', 'Centimeters');
    pos = get(pdf_plot, 'Position');
    set(pdf_plot, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Centimeters', ...
        'PaperSize', [pos(3), pos(4)]);
    % print(pdf_plot, [filename, '.pdf'], '-dpdf', '-r0'); % Full quality
    print(pdf_plot, [filename, '.pdf'], '-dpdf', ['-r', num2str(dpi)]);

    % Save png
    % print(pdf_plot, [filename, '.png'], '-dpng', '-r600'); % Good quality
    print(pdf_plot, [filename, '.png'], '-dpng', ['-r', num2str(dpi)]);

end
