function main( verbosity )
%
% Here we collect the numerical experiments run for the bamphi manuscript.
% For any question contact:
%              Franco Zivcovich, franco.zivcovich@gmail.com
%
% if verbosity = 0, we display almost nothing
% if verbosity > 0, we display some information
% if verbosity = 1, we also display a short description of the experiment
%
  if not( nargin )
    verbosity = 1;
  end
  clc
  if verbosity == 1
    fprintf('If you want me to talk less, set verbosity lower. Valid verbosity options are:\n')
    fprintf('    - 0 for shut up;\n')
    fprintf('    - anything between 0 and 1 for few comments;\n')
    fprintf('    - 1 for me talking a lot.\n\n')
  end
  addpath( genpath('../') )

  % Choose the tests you want to run
  % (just comment the lines of the tests you DON'T want to run):
  which_tests{ 1 } = 'adr'; % Advection-Diffusion-Reaction equation
  which_tests{ 2 } = 'brg'; % Burger equation
  which_tests{ 3 } = 'sch_S'; % smooth cubic Schrodinger equation
  which_tests{ 4 } = 'sch_L'; % low-regularity cubic Schrodinger equation
  which_tests{ 5 } = 'skg'; % low-regularity cubic Schrodinger equation
  which_tests{ 6 } = 'st_ac'; % Stochastic Allen-Cahn


  for test_idx = 1 : length( which_tests )

    if strcmp( which_tests{ test_idx }, 'adr' )
      label = 'adr';
      test.equation      = 'adr';
      test.Ns            = [ 500 500 ];
      test.Nt            = 2^8 : 2^8 : 2^11;
      test.integrator    = 'exprk4s6';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end

    if strcmp( which_tests{ test_idx }, 'brg' )
      label = 'brg';
      test.equation      = 'brg';
      test.Ns            = [ 500, 500 ];
      Nt = 2.^[4:7];
      Nt = sort( [Nt, cumsum( Nt ) / 2 ]);
      test.Nt            = Nt( 1:end-2 );
      test.integrator    = 'epirk4s3a';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end

    if strcmp( which_tests{ test_idx }, 'sch_S' )
      label = 'sch_S';
      test.equation      = 'sch_S';
      test.Ns            = [ 200, 200 ];
      test.Nt            = 2.^[ 1 : 9 ];
      test.integrator    = 'schss68';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end

    if strcmp( which_tests{ test_idx }, 'sch_L' )
      label = 'sch_L';
      test.equation      = 'sch_L';
      test.Ns            = [ 200 200 ];
      test.Nt            = 2.^[ 2 : 11 ];
      test.theta         = 5 / 2;
      test.integrator    = 'schexpli2s4';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end

    if strcmp( which_tests{ test_idx }, 'skg' )
      % This test is messy 'cause we launch kiops in a way and bamphi in another.
      % MODALITIES (for more details, see BAMPHI paper):
      % AA: I use the (very) sparse, highly non-normal A to do everything
      % RR: I use the  less  sparse, perfectly  normal R to do everything
      % RA: I compute Ritz's values with R and run calculations with A (only bamphi)
      label = 'skg';
      test.equation      = 'skg';
      test.Ns            = [ 500 500 ];
      test.Nt            = 2.^[ 4 : 11 ];
      test.integrator    = 'exptg2s1';
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %

      % kiops
      test.routines      = {'kiops'};
      test.mode = 'AA';
      memo = feval( test.equation, test );
      store.(label).memo.(test.integrator).(test.routines{1}) = memo.(test.integrator).(test.routines{1});
      kiops_testmode = [ 'kiops', test.mode ];

      % bamphi
      test.routines      = {'bamphi'};
      test.mode = 'RA';
      memo = feval( test.equation, test );
      store.(label).memo.(test.integrator).(test.routines{1}) = memo.(test.integrator).(test.routines{1});
      bamphi_testmode = [ 'bamphi', test.mode ];

      % wrap up
      test.routines = {'bamphi','kiops'};
      test.mode     = {kiops_testmode,bamphi_testmode};
      store.(label).test = test;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end


    if strcmp( which_tests{ test_idx }, 'st_ac' )
      label = 'st_ac';
      test.equation      = 'st_ac';
      test.Ns            = [ 200 200 ];
      test.Nt            = 2.^[ 1 : 8 ]; % important: timesteps must be one multiple of the other for stochasticity function
      test.MoCa          = 100;%1e2;
      test.integrator    = 'setdm1';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      %
      if verbosity == 1
        help (which_tests{ test_idx })
        fprintf('\n')
        help ([test.integrator,'_bamphi'])
        fprintf('\n')
        help ([test.integrator,'_kiops'])
        fprintf('\n DISCLAIMER: \n If it''s taking too long you can modify the space and time discretization by \n changing test.Ns and test.Nt in main.m.\n\n')
      end
      %
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test store
    end
    if verbosity == 1
      fprintf('\n\n\n\n\n')
    else
      fprintf('\n')
    end
  end
   rmpath( genpath('../') )
end

function final_display( test )

  fprintf(['\n'])
  fprintf(['To access the results use load(''./results/',test.equation,'.mat''), you''ll find\n'])
  fprintf(['    the test information calling store.',test.equation,'.test\n'])
  fprintf(['    the test     results calling store.',test.equation,'.memo.',test.integrator,'.bamphi \n'])
  fprintf(['                              or store.',test.equation,'.memo.',test.integrator,'.kiops \n'])
  fprintf(['Also, you can plot results by calling\n'])
  fprintf(['    test_plotter(''',test.equation,''',''precision'' ) for the precision plot\n'])
  fprintf(['    test_plotter(''',test.equation,''',''order''     ) for the order     plot\n'])
  fprintf(['                                                                       :)\n'])

end
