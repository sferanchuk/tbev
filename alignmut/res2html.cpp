
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main( int argc, char **argv )
{
	const int nkeys = 5;
	const char *key[] = { "<Pos>", "<SPos>", "<Letter>", "<Percent>", "<Other>" };
	char buf[4096];
	printf( "<html><body><table>\n" );
	while ( fgets( buf, 4096, stdin ) )
	{
		if ( !strstr( buf, "<Line>" ) ) continue;
		printf( "<tr>" );
		for ( int kc = 0; kc < nkeys; kc++ )
		{
			char *v = strstr( buf, key[kc] ) + strlen( key[kc] );
			printf( "<td>%.*s</td>", strcspn( v, "<" ), v );
		}
		printf( "</tr>\n" );
	}
	printf( "</table></body></html>\n" );
}
		
