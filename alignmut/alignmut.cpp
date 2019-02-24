
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;

int main( int argc, char **argv )
{
	if ( argc == 1 )
	{
		fprintf( stderr, "arguments: beg end [size(20)]\n" );
		return 1;
	}
	const char *anames[] = { "TBEV", "OHFV", "KFDV", "YFV", "JEV", "WNV", "DV1", "DV2", "DV3", "DV4" };
	int off = 0;
	int beg = atoi( argv[1] ) - off;
	int end = atoi( argv[2] ) - off;
	int osize = 20;
	if ( argc > 3 ) osize = atoi( argv[3] ); 
	vector<map<int,pair<string,string> > > lmap;
	for ( int i = 1; i <= 10; i++ )
	{
		char nbuf[256];
		sprintf( nbuf, "res_%d.xml", i );
		FILE *ifile = fopen( nbuf, "rt" );
		if ( !ifile ) break;
		lmap.resize( lmap.size() + 1 );
		char buf[4096];
		while ( fgets( buf, 4096, ifile ) )
		{
			char *p1 = strstr( buf, "<SPos>" );
			if ( !p1 ) continue;
			string s1( p1 + 6, 0, strcspn( p1 + 6, "<" ) );
			char *p2 = strstr( buf, "<Letter>" );
			string s2( p2 + 8, 0, strcspn( p2 + 8, "<" ) );
			char *p3 = strstr( buf, "<Percent>" );
			string s3( p3 + 9, 0, strcspn( p3 + 9, "<" ) );
			char *p4 = strstr( buf, "<Other>" );
			string s4( p4 + 7, 0, strcspn( p4 + 7, "<" ) );
			char *p5 = strstr( buf, "<LCount>" );
			string s5( p5 + 8, 0, strcspn( p5 + 8, "<" ) );
			int pos = atoi( s1.data() );
			int LCount = atoi( s5.data() );
			if ( pos < beg || pos > end ) continue;
			multimap<double,string> let;
			char lbuf[256];
			strcpy( lbuf, s4.data() );
			double xp = 0;
			for ( char *p = strtok( lbuf, " " ); p; p = strtok( 0, " \t" ) )
			{
				string clet = p;
				p = strtok( 0, " \t" );
				let.insert( pair<double,string>( atof( p ), clet ) );
				if ( clet == "X" ) xp += atof( p );
			}
			let.insert( pair<double,string>( atof( s3.data() ), s2 ) );
			if ( s2 == "X" ) xp += atof( s3.data() );
			string l1 = let.rbegin()->second;
			if ( ( let.rbegin()->first + xp ) * LCount * 0.01 < LCount - 2.2 ) l1 = "<font color=\"green\">" + let.rbegin()->second + "</font>";
			pair<string,string> cp( l1, "&nbsp;" );
			double entropy = 0;
			for ( multimap<double,string>::iterator it = let.begin(); it != let.end(); it++ ) entropy += 0.01 * it->first * log( 0.01 * it->first );
			multimap<double,string>::reverse_iterator it = let.rbegin();
			it++;
			if ( it != let.rend() ) 
			{
				cp.second = it->second;
				if ( entropy < -0.5 && let.size() > 2 ) cp.second = "<font color=green>#</font>";
				else if ( it->first * LCount * 0.01 > 1.5 ) cp.second = "<font color=green>" + it->second + "</font>";
			}
			lmap.back()[ pos ] = cp;
		}
		fclose( ifile );
	}
	printf( "<html><body bgcolor=white><font face=\"courier\" size=\"10px\">\n" );
	for ( int k0 = beg; k0 <= end; k0 += osize )
	{
		printf( "<table style=\"{ margin-top:0px; margin-bottom:0px; }\"><tr><td>\n" );
		char strnum[2][16];
		sprintf( strnum[0], "%d", k0 + off );
		sprintf( strnum[1], "%d", min( end, k0 + osize ) + off );
		string title( strnum[0] );
		for ( int i = 0; i < min( k0 + osize, end ) - k0 - 1; i++ )
		{
			//if ( ( i - beg  ) % 10 == 0 ) title.append( "|<sub>&nbsp;</sub><font size=\"1px\">&nbsp;</font>" );
			//else 
			title.append( "&nbsp;<sub>&nbsp;</sub><font size=\"1px\">&nbsp;</font>" );
		}
		title.append( strnum[1] );
		printf( "<td>%s", title.data() );
		printf( "\n" );
		for ( int i = 0; i < lmap.size(); i++ )
		{
			printf( "<tr><td><b>%s&nbsp;&nbsp;</b><td><b>", anames[i] );
			int k = 0;
			for ( map<int,pair<string,string> >::iterator it = lmap[i].begin(); it != lmap[i].end(); it++ )
			{
				if ( it->first < k0 || it->first >= k0 + osize ) continue;
				printf( "%s<font color=\"red\"><sub>%s</sub></font><font size=\"1px\">&nbsp;</font>", it->second.first.data(), it->second.second.data() );
				k++;
				//if ( k % 10 == 0 ) printf( "&nbsp;<td><b>" );
			}
			printf( "</b>\n" );
		}
		printf( "</table></b>" );//<hr style=\"height:5pt; visibility:hidden;\" />" );
	}
	printf( "</font></body></html>\n" );
	return 0;
}