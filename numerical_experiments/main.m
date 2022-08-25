function main( verbosity )
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
  % (just comment the lines with the tests you DON'T want to run):
  which_tests{ 1 } = 'adr'; % Advection-Diffusion-Reaction equation
  which_tests{ 2 } = 'brg'; % Burger equation
  which_tests{ 3 } = 'scS'; % smooth cubic Schrodinger equation
  which_tests{ 4 } = 'sch'; % low-regularity cubic Schrodinger equation
  which_tests{ 5 } = 'skg'; % low-regularity cubic Schrodinger equation
  which_tests{ 6 } = 'sac'; % Stochastic Allen-Cahn


  for test_idx = 1 : length( which_tests )

    if strcmp( which_tests{ test_idx }, 'adr' )
      if verbosity == 1
        help adr
      end
      label = 'adr';
      test.equation      = 'adr';
      test.Ns            = [ 500 500 ] / 50;
      test.Nt            = 2^8 : 2^8 : 2^11;
      test.integrator    = 'exprk4s6';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test
    end

    if strcmp( which_tests{ test_idx }, 'brg' )
      if verbosity == 1
        help brg
      end
      label = 'brg';
      test.equation      = 'brg';
      test.Ns            = [ 500, 500 ] / 50;
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
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test
    end

    if strcmp( which_tests{ test_idx }, 'scS' )
      if verbosity == 1
        help scS
      end
      label = 'scS';
      test.equation      = 'scS';
      test.Ns            = [ 100, 100 ] / 10;
      test.Nt            = 2.^[ 1 : 8 ];
      test.integrator    = 'cuschss68';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test
    end

    if strcmp( which_tests{ test_idx }, 'sch' )
      if verbosity == 1
        help sch
      end
      label = 'sch';
      test.equation      = 'sch';
      test.Ns            = [ 100 100 ] / 10;
      test.Nt            = 2.^[ 3 : 11 ];
      test.theta         = 5 / 2;
      test.integrator    = 'cuschRS21';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test
    end

    if strcmp( which_tests{ test_idx }, 'skg' )
      if verbosity == 1
        help skg
      end
      % This test is messy 'cause we launch kiops in a way and bamphi in another.
      % MODALITIES (for more details, see BAMPHI paper):
      % AA: I use the (very) sparse, highly non-normal A to do everything
      % RR: I use the  less  sparse, perfectly  normal R to do everything
      % RA: I compute Ritz's values with R and run calculations with A (only bamphi)
      label = 'skg';
      test.equation      = 'skg';
      test.Ns            = [ 500 500 ] / 50;
      test.Nt            = 2.^[ 4 : 11 ];
      test.integrator    = 'exptg2s1';
      test.compute_error = false;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end

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

      clear memo label test
    end


    if strcmp( which_tests{ test_idx }, 'sac' )
      if verbosity == 1
        help sac
      end
      label = 'sac';
      test.equation      = 'sac';
      test.Ns            = [ 100 100 ] / 10;
      test.Nt            = 2.^[ 2 : 11 ]; % important: timesteps must be one multiple of the other for stochasticity function
      test.MoCa          = 1e2;
      test.integrator    = 'stochexprk1';
      test.routines      = {'bamphi','kiops'};
      test.compute_error = true;
      %
      if verbosity > 0
        fprintf('Test launched on '),fprintf('\b'),disp(datetime(now,'ConvertFrom','datenum'))
        disp('Test information:'),disp(test)
      end
      memo = feval( test.equation, test );
      store.(label).test = test;
      store.(label).memo = memo;
      save( [ './results/', label ], 'store' )

      if verbosity > 0
        final_display( test )
      end

      clear memo label test
    end
    fprintf('\n\n\n\n\n')
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
