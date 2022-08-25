function [ D1, D2, h, varargout ] = space_discr( omega, N, boundary )
% Returns the D1 and D2 Finite Differences matrices, the spacings h and, if
% requested, the d-dimensional ndgrid of the dimensions' linspaces.
% The inputs are the domain omega, the number N of discretization points for each
% dimension and the kind of boundary conditions to set.

  space_vars = fieldnames( omega );
  d = length( N );
  for i = 1 : d
    if strcmp( boundary.(space_vars{ i }).l, 'hom_dirichlet' )
      if strcmp( boundary.(space_vars{ i }).r, 'hom_dirichlet' )
        lin{ i } = linspace( omega.(space_vars{ i }).l, omega.(space_vars{ i }).r, N( i ) + 2 );
        lin{ i } = lin{ i }( 2 : end - 1 );
        h.(space_vars{i}) = ( omega.(space_vars{ i }).r - omega.(space_vars{ i }).l ) / ( N( i ) +  1 );
      elseif strcmp( boundary.(space_vars{ i }).r, 'hom_neumann' )
        lin{ i } = linspace( omega.(space_vars{ i }).l, omega.(space_vars{ i }).r, N( i ) + 1 );
        lin{ i } = lin{ i }( 2 : end );
        h.(space_vars{i}) = ( omega.(space_vars{ i }).r - omega.(space_vars{ i }).l ) / ( N( i ) + 0 );
      else
        error('Somethings wrong :C');
      end
    elseif strcmp( boundary.(space_vars{ i }).l, 'hom_neumann' )
      if strcmp( boundary.(space_vars{ i }).r, 'hom_neumann' )
        lin{ i } = linspace( omega.(space_vars{ i }).l, omega.(space_vars{ i }).r, N( i ) );
        lin{ i } = lin{ i }( 1 : end );
        h.(space_vars{i}) = ( omega.(space_vars{ i }).r - omega.(space_vars{ i }).l ) / ( N( i )  - 1 );
      elseif strcmp( boundary.(space_vars{ i }).r, 'hom_dirichlet' )
        lin{ i } = linspace( omega.(space_vars{ i }).l, omega.(space_vars{ i }).r, N( i ) + 1 );
        lin{ i } = lin{ i }( 1 : end - 1 );
        h.(space_vars{i}) = ( omega.(space_vars{ i }).r - omega.(space_vars{ i }).l ) / ( N( i ) + 0 );
      else
        error('Somethings wrong :C');
      end
    else
      error('Somethings wrong :C');
    end
    % D1
    v0 = ones( N( i ),1 ) / ( 2 * h.(space_vars{i}) );
    D1.(space_vars{ i }) = spdiags( v0 * [ -1,1 ],[ -1, 1 ], N( i ), N( i ) );
    if strcmp( boundary.(space_vars{ i }).l, 'hom_neumann' )
      D1.(space_vars{ i })( 1, 2 ) = 0;
    end
    if strcmp( boundary.(space_vars{ i }).r, 'hom_neumann' )
      D1.(space_vars{ i })( N( i ), N( i ) - 1 ) = 0;
    end
    % D2
    v0 = ones( N( i ),1 ) / ( h.(space_vars{i})^2 );
    D2.(space_vars{ i }) = spdiags( v0 * [ 1, -2, 1 ], [ -1, 0, 1 ], N( i ), N( i ) );
    if strcmp( boundary.(space_vars{ i }).l, 'hom_neumann' )
      D2.(space_vars{ i })( 1, 2 ) = 2 / ( h.(space_vars{i})^2 );
    end
    if strcmp( boundary.(space_vars{ i }).r, 'hom_neumann' )
      D2.(space_vars{ i })( N( i ), N( i ) - 1 ) = 2 / ( h.(space_vars{i})^2 );
    end
  end
  varargout = cell( d,1 );
  [ varargout{:} ] = ndgrid( lin{:} );
end
