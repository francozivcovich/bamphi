function test_plotter( test2plot, plot_kind )

  % settings
  if nargin < 2
    plot_kind = 'precision';
  end
  marker.bamphi = '-o';
  marker.kiops  = '-x';
  save_fig  = true; path2png   = './results/png-figs/';
  % save_tikz = true; path2tikz  = './results/tikz-figs/';
  % if save_tikz % Sorry: atm I don't have time to find a different way to do it.
  %   addpath('/Users/francozeta/work/matlab2tikz-master/src') % just comment this and turn save_tikz to false pls
  % end

  % recover data
  filePattern = fullfile( './results/', '*.mat' );
  matFiles = dir( filePattern );
  for i = 1 : length( matFiles )
    dataName{ i } = matFiles( i ).name;
    load( [ './results/', dataName{ i } ] )
    test = fieldnames( store );
    for j = 1 : length( test )
      if not( nargin ) || isempty( test2plot )
        superstore.(test{ j }) =  store.(test{ j });
      else
        if strcmp( test2plot, test{ j } )
          superstore.(test{ j }) =  store.(test{ j });
        end
      end
    end
  end

  % plot
  equation = fieldnames( superstore );
  for e = 1 : length( equation )
    testinfo = superstore.(equation{e}).test;
    testmemo = superstore.(equation{e}).memo;
    integrator = fieldnames( testmemo );
    for i = 1 : length( integrator )
      routine = fieldnames( testmemo.(integrator{i}) );
      idx = 1;
      figure,
      for r = 1 : length( routine )
        error = testmemo.(integrator{i}).(routine{r}).error;
        matve = testmemo.(integrator{i}).(routine{r}).matve;
        clock = testmemo.(integrator{i}).(routine{r}).clock;
        tstep = testmemo.(integrator{i}).(routine{r}).tstep;

        if strcmp( plot_kind, 'precision' )
          semilogy( clock, error, marker.(routine{r}) ), hold on
        end
        if strcmp( plot_kind, 'order' )
          loglog( tstep, error, marker.(routine{r}) ), hold on
        end
        leg{ idx } = [ routine{r}, ' ', integrator{i} ]; idx = idx + 1;
      end
    end
    legend( leg{:},'Location','NorthEast' )
    % title and filename
    my_title    = [ equation{e},' ',integrator{i} ];
    my_filename = [ equation{e},    integrator{i} ];
    title( my_title )
    if save_fig
      saveas( gcf, [ path2png, my_filename,'_', plot_kind, '.png' ] ),
    end
    % if save_tikz
    %   matlab2tikz( [ path2tikz, my_filename,'_', plot_kind, '.tikz' ], 'showInfo', false ),
    % end
    clear leg
    hold off
  end

  % if save_tikz
  %   rmpath('/Users/francozeta/work/matlab2tikz-master/src')
  % end

end
