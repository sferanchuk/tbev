
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main( int argc, char **argv )
{
	FILE *ifile = stdin;
	int vc = 0;
	double sum = 0;
	double x[3] = { 0 };
	char buf[ 256 ];
	
	while ( fgets( buf, 256, ifile ) )
	{
		const char *sep = "c(, \t\n";
		for ( char *p = strtok( buf, sep ); p; p = strtok( 0, sep ) )
		{
			double v = atof( p );
			int n1 = vc % 3;
			int n2 = ( vc / 3 ) % 10;
			int n3 = vc / 30;
			if ( n1 == 0 )
			{
				if ( n2 == 0 && n3 > 0 )
				{
					printf( "s_%d_%d\t%g\t\n", ( n3 - 1 ) * 10 + 1, n3 * 10 + 1, 1./ sum );
					sum = 0;
				}
				sum += sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
				x[0] = x[1] = x[2] = 0;
			}
			x[ n1 ] = v;
			vc++;
			//fprintf( stderr, "%d %d %d %d\n", vc, n1, n2, n3 );
		}
	}
	return 0;
}

	