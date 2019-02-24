
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>
#include <vector>
#include <set>
#include <map>

using namespace std;

int main( int argc, char **argv )
{
	vector< vector< string > > table;
	vector<int> mut;
	vector< string > ids;
	vector< map<string,int> > other;
	vector<int> lcvec;
	int fcnt = 1;
	do
	{
		char nbuf[32];
		sprintf( nbuf, "res_%d.xml", fcnt );
		fcnt++;
		FILE *ifile = fopen( nbuf, "rt" );
		if ( !ifile ) break;
		const int nkeys = 7;
		const char *key[] = { "<Pos>", "<SPos>", "<Letter>", "<Percent>", "<LCount>", "<Other>", "<Ids>" };
		char buf[4096];
		int lc = 0;
		while ( fgets( buf, 4096, ifile ) )
		{
			if ( !strstr( buf, "<Line>" ) ) continue;
			vector<string> fields;
			for ( int kc = 0; kc < nkeys; kc++ )
			{
				char *v = strstr( buf, key[kc] ) + strlen( key[kc] );
				fields.push_back( string( v, 0, strcspn( v, "<" ) ) );
			}
			if ( table.size() <= lc ) 
			{
				table.resize( lc + 1 );
				ids.resize( lc + 1 );
				other.resize( lc + 1 );
				lcvec.resize( lc + 1, 0 );
				table[lc].push_back( fields[0] );
				table[lc].push_back( fields[1] );
				table[lc].push_back( fields[2] );
			}
			table[lc].push_back( fields[3] );
			if ( mut.size() <= lc ) mut.resize( lc + 1, 0 );
			int lcount = atoi( fields[4].data() );
			if ( lcount == 0 ) mut[ lc ] = 10;
			else mut[lc] += int(  lcount * ( 100 - atof( fields[3].data() ) ) * 0.01 + 0.3 );
			ids[lc].append( fields[6] );
			char nbuf[256];
			strcpy( nbuf, fields[5].data() );
			for ( char *p = strtok( nbuf, " \t" ); p; p = strtok( 0, " \t" ) )
			{
				string clet = p;
				double f = atof( strtok( 0, " \t" ) );
				other[lc][ clet ] += f * lcount * 0.01 + 0.4;
			}
			lcvec[lc] += lcount;
			lc++;
		}
	}
	while ( 1 );
	printf( "<html><body><table>" );
	for ( int lc = 0; lc < table.size(); lc++ )
	{
		printf( "<tr>" );
		for ( int fc = 0; fc < table[lc].size(); fc++ ) printf( "<td>%s", table[lc][fc].data() );
		char symb[4] = ".";
		if ( mut[lc] == 0 ) strcpy( symb, "*" );
		else if ( mut[lc] < 3 ) sprintf( symb, "%d", mut[lc] );
		printf( "<td>%s<td>", symb );
		multimap<int,string> rother;
		for ( map<string,int>::iterator it = other[lc].begin(); it != other[lc].end(); it++ ) rother.insert( pair<int,string>( it->second, it->first ) );
		for ( multimap<int,string>::reverse_iterator it = rother.rbegin(); it != rother.rend(); it++ ) printf( "%s %d; ", it->second.data(), it->first );
		printf( "<td>" );
		set<string> cids;
		char buf[4096];
		strcpy( buf, ids[lc].data() );
		set<string> sids;
		for ( char *p = strtok( buf, " " ); p; p = strtok( 0, " " ) ) sids.insert( p );
		int ocnt = 0;
		for ( set<string>::iterator it = sids.begin(); it != sids.end() && ocnt < 3; it++, ocnt++ ) 
		{
			vector<string> idi;
			char nbuf[256];
			strcpy( nbuf, it->data() );
			for ( char *p = strtok( nbuf, "|" ); p; p = strtok( 0, "|" ) ) idi.push_back( p );
			printf( "%s ", idi.back().data() );
		}
		printf( "\n" );
	}
	printf( "</table></body></html>" );
	return 0;		
}
		
