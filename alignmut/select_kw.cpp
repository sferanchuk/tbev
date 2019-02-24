
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main( int argc, char **argv )
{
    const char *kw = argv[2];
    FILE *ifile = fopen( argv[1], "rt" );
    char buf[65535];
    while ( fgets( buf, 65535, ifile ) )
    {
	if ( buf[0] != '>' ) continue;
	if ( !strstr( buf, kw ) ) continue;
	printf( "%s", buf + 1 );
    }
    fclose( ifile );
    return 0;
}