
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <string>
#include <map>
#include <list>
#include <set>

using namespace std;

#include <pdb++.h>

#include "vector.h"
#include "geom-f.h"
#include "format.h"
#include "indicator.h"
#include "simrec.h"
#include "cluster.h"
#include "motifs.h"
#include "motifs3d.h"
#include "distmatr.h"

#include "rmsdstats.cpp"

const double MaxDist = 1.8;
const double MinPower = 3;
static vector<double> cMinPower;
static int MaxNumClusters = 15;//8;
static int MaxNumPClusters = 50;
static int MaxNumInCluster = 4;
const int MaxNumSS = 15;

extern bool connect_pair_clusters; // from wmotif.cpp
SCluster *Merge( const SCluster *m1, const SCluster *m2, int width );
void BuildBeta( VBB& chain, int length );
void BuildAlpha( VBB& chain, int length );

typedef pair<int,int> IPair;

struct StatLetter
{
	int pos;
	char letter;
	double distr[26];
};

typedef vector<StatLetter> StatSequence;

struct TemplateStr
{
	string sequence;
	string ss;
	VBB bb;
	vector<NAtom> natom;
	vector<SCluster> cluster;
	int order;
	int type;
	int chain;
	int id;
};

static const char *aaLetters = "ACDEFGHIKLMNPQRSTVWY";

static double aa_frequences[ 26 ] = {
	0.085, 0, 0.018, 0.055, 0.066, 0.040,              // A,B,C,D,E,F
	0.075, 0.022, 0.059, 0, 0.059, 0.094,              // G,H,I,J,K,L
	0.025, 0.044, 0, 0.052, 0.040, 0.050,              // M,N,P,Q,R
	0.069, 0.060, 0, 0.072, 0.014, 0, 0.034, 0          // S,T,V,W,Y,-
	};

#include "pam30.cpp"

double mmean;
double mmdisp;

static void NormalizeWM( double& mean, double& mdisp )
{
	mean = 0;
    for ( int i = 0; i < 26; i++ ) 
	{
		for ( int j = 0; j < 26; j++ ) 
		{
			mean += aa_frequences[i] * aa_frequences[j] * weightMatrix[i][j];
		}
	}
	double dispersion = 0;
    for ( int i = 0; i < 26; i++ ) 
	{
		for ( int j = 0; j < 26; j++ ) 
		{
			dispersion += aa_frequences[i] * aa_frequences[j]
				* ( weightMatrix[i][j] - mean )
				* ( weightMatrix[i][j] - mean );
		}
	}
	mdisp = sqrt( dispersion );
}


int clusterize( int size, double *_distances, double threshold, void *pStr, int *sizes, int **members );

static double CalcMeasure( double rmsd, int size )
{
	if ( size > MaxSegmLength ) size = MaxSegmLength;
	double power = RmsdParams[ size ][ 0 ];
	double lambda = RmsdParams[ size ][ 1 ];
	double x = rmsd * lambda;
	if ( x > power ) return 0;
	double rv = -log( x ) * power + x + gamma( power + 1 );
	return rv;
}

void ProcessSeqDothelix( StatSequence& seq, TemplateStr& ts )
{
	vector<int> begs;
	begs.push_back( 0 );
	for ( int bc = 1; bc < ts.bb.size() - 1; bc++ )
	{
		if ( ( ts.bb[ bc ].v[ aC ] - ts.bb[ bc + 1 ].v[ aN ] ).norm() > MaxDist ) begs.push_back( bc );
	}
	string& cseq = ts.sequence;
	double PowerThreshold = 3.5;
	vector<Homology> hom;
	for ( int sc = -( (int)( seq.size() ) - 1 ); sc < int( cseq.size() ); sc++ )
	{
		int beg1 = max( 0, -sc );
		int beg2 = max( 0, sc );
		vector<double> elem;
		for ( int lc = 0; beg1 + lc < seq.size() && beg2 + lc < cseq.size(); lc++ )
		{
			char cletter = cseq[ beg2 + lc ];
			double v;
			if ( cletter == 'U' || cletter == '-' ) v = -20;
			else v = seq[ beg1 + lc ].distr[ cletter - 'A' ];
			elem.push_back( v );
		}
		for ( int bc = 0; bc < begs.size(); bc++ )
		{
			int cbeg2 = max( begs[bc], beg2 );
			int cend2 = ( bc + 1 == begs.size() ) ? int( elem.size() ) + beg2 - cbeg2 : min( max( begs[ bc + 1 ] - cbeg2, 0 ), int( elem.size() ) + beg2 - cbeg2 );
			if ( cend2 <= 0 ) continue;
			//if ( cend2 <= cbeg2 ) continue;
			ushort num_found = 0;
			Homology_info *hi = DotHelix( &(elem[cbeg2 - beg2]), cend2, 4, MinPower, num_found );
			for ( ushort hc = 0; hc < num_found; hc++ )
			{
				Homology h;
				h.beg.resize(2);
				h.x = h.beg[0] = beg1 + cbeg2 - beg2 + hi[hc].left;
				h.y = h.beg[1] = cbeg2 + hi[hc].left;
				h.length = hi[hc].length;
				h.sum = hi[hc].sum;
				h.next = 0;
				hom.push_back( h );
			}
			//if ( num_found > 0 ) printf( "%d %d %d %d %d %d %d %d\n", ts.order, bc, begs[bc], beg2, elem.size(), cbeg2, cend2, num_found );
			delete hi;
		}
	}
	//FindClusters( hom, ts.cluster, 3 );
	if ( hom.size() )
	{
		ts.cluster.resize( 1 );
		FindBestWay( hom, seq.size(), cseq.size(), ts.cluster[0] );
		if ( ts.cluster[0].power != ts.cluster[0].power ) ts.cluster.erase( ts.cluster.begin() );
	}
	for ( int oc = 0; oc < ts.cluster.size(); oc++ )
	{
		ts.cluster[oc].order.resize( 2 );
		ts.cluster[oc].order[0] = ts.chain;
		ts.cluster[oc].order[1] = ts.order;
		for ( int hc = 0; hc < ts.cluster[oc].h.size(); hc++ )
		{
			ts.cluster[oc].h[hc].beg.resize( 2 );
			ts.cluster[oc].h[hc].beg[0] = ts.cluster[oc].h[hc].x;
			ts.cluster[oc].h[hc].beg[1] = ts.cluster[oc].h[hc].y;
		}
	}
	/*
	double pmax = 0;
	for ( int cc = 0; cc < ccluster.size(); cc++ )
	{
		if ( ccluster[cc].power > pmax )
		{
			pmax = ccluster[cc].power;
			cluster[0] = ccluster[cc];
		}
	}
	*/
	//if ( hom.size() ) FindBestWay( hom, seq.size(), cseq.size(), cluster[0] );
}

double CompareTemplates( TemplateStr& t1, TemplateStr& t2, SClusterSet& scset )
{
	double mmax = 0;
	for ( int cc1 = 0; cc1 < t1.cluster.size(); cc1++ )
	{
		for ( int cc2 = 0; cc2 < t2.cluster.size(); cc2++ )
		{
			SCluster& cl1 = t1.cluster[cc1];
			SCluster& cl2 = t2.cluster[cc2];
			if ( cl1.power < cMinPower[cl1.order[0]] || cl2.power < cMinPower[cl2.order[0]] ) continue;
			SCluster clr;
			clr.order.resize( 3 );
			clr.order[0] = t1.chain;
			clr.order[2] = t1.order;
			clr.order[1] = t2.order;
			clr.sum = 0;
			vector<Vector> v1;
			vector<Vector> v2;
			for ( int hc1 = 0; hc1 < cl1.h.size(); hc1++ )
			{
				Homology& h1 = cl1.h[hc1];
				for ( int hc2 = 0; hc2 < cl2.h.size(); hc2++ )
				{
					Homology& h2 = cl2.h[hc2];
					Homology hr;
					hr.beg.resize( 3 );
					if ( h2.x + h2.length < h1.x ) continue;
					if ( h2.x > h1.x + h1.length ) continue;
					hr.length = 0;
					hr.sum = 0;
					for ( int pc = max( h1.x, h2.x ); pc < min( h1.x + h1.length, h2.x + h2.length ); pc++ )
					{
						if ( hr.length == 0 )
						{
							hr.beg[0] = pc; 
							hr.beg[2] = h1.y + pc - h1.x;
							hr.beg[1] = h2.y + pc - h2.x;
						}
						v1.push_back( t1.bb[ h1.y + pc - h1.x ].v[1] );
						v2.push_back( t2.bb[ h2.y + pc - h2.x ].v[1] );
						hr.length++;
						hr.sum += ( weightMatrix[ t1.sequence[ h1.y + pc - h1.x ] - 'A' ][ t2.sequence[ h2.y + pc - h2.x ] - 'A' ] - mmean ) / mmdisp;
					}
					if ( hr.length > 0 )
					{
						hr.sum += ( h1.sum * hr.length ) / h1.length + ( h2.sum * hr.length ) / h1.length;
						clr.h.push_back( hr );
						clr.sum += hr.sum;
						clr.length += hr.length;
					}
				}
			}
			if ( v1.size() < 5 ) continue;
			clr.length = clr.h.back().beg[0] + clr.h.back().length - clr.h[0].beg[0];
			double rmsd = FindRmsd( v1.size(), &(v1[0]), &(v2[0]) );
			clr.power = clr.sum / sqrt( 3 * clr.length );
			//if ( rmsd > 5 ) continue;
			double measure = CalcMeasure( rmsd, v1.size() ) * sqrt( double( v1.size() ) / clr.length );
			double coeff = 1;
			/* = ( fmin( 
				-log( 1. - ( 1. - exp( -fmin(measure,20.) ) ) * ( 1.-exp( -fmin(clr.power,20)) ) ),
				fmin( cl1.power, cl2.power ) ) ) / clr.power;
			*/
			for ( int hc = 0; hc < clr.h.size(); hc++ ) 
			{
				clr.h[hc].sum *= coeff;
			}
			clr.sum *= coeff;
			clr.power *= coeff;
			if ( clr.power > 100 )
			{
				printf( "karaul\n" );
			}
			if ( clr.power >= MinPower && measure > 2 )	
			{
				SCluster *tcl1 = Merge( &clr, &cl1, clr.order.back() + 1 );
				if ( tcl1 )
				{
					SCluster *tcl2 = Merge( tcl1, &cl2, clr.order.back() + 1 );
					delete tcl1;
					if ( tcl2 )
					{
						scset.cluster.push_back( *tcl2 );
						delete tcl2;
					}
				}
			}
		}
	}
	return mmax;
}

bool CompareTemplates2( TemplateStr& t1, TemplateStr& t2, SCluster& cl, displacement& displ )
{
	int o1 = -1;
	int o2 = -1;
	for ( int oc = 0; oc < cl.order.size(); oc++ )
	{
		if ( cl.order[oc] == t1.order ) o1 = oc;
		if ( cl.order[oc] == t2.order ) o2 = oc;
	}
	if ( o1 == -1 || o2 == -1 ) return false;
	vector<Vector> v1;
	vector<Vector> v2;
	for ( int hc = 0; hc < cl.h.size(); hc++ )
	{
		Homology& h = cl.h[hc];
		if ( h.beg[o1] == -1 || h.beg[o2] == -1 ) continue;
		for ( int pc = 0; pc < h.length; pc++ )
		{
			v1.push_back( t1.bb[ h.beg[o1] + pc ].v[1] );
			v2.push_back( t2.bb[ h.beg[o2] + pc ].v[1] );
		}
	}
	if ( v1.size() == 0 ) return false;
	double rmsd;
	displacement *d0 = FindDisplacement( v1.size(), &(v1[0]), &(v2[0]), &rmsd );
	if ( !d0 ) return false;
	displ = *d0;
	delete d0;
	return true;
}

void DoTransform( TemplateStr& t, displacement& d, VBB& ncoord, vector<NAtom>& natom )
{
	for ( int rc = 0; rc < t.bb.size(); rc++ )
	{
		for ( int ac = 0; ac < 4; ac++ )
		{
			ncoord[rc].v[ac] = d.Move( t.bb[rc].v[ac] );
		}
	}
	for ( int ac = 0; ac < t.natom.size(); ac++ )
	{
		NAtom na = t.natom[ ac ];
		na.coord = d.Move( na.coord );
		natom.push_back( na );
	}
}

void WriteDisplacement( const displacement& d, Line& l )
{
	Vector vx( d.Turn.xx, d.Turn.xy, d.Turn.xz );
	Vector vy( d.Turn.yx, d.Turn.yy, d.Turn.yz );
	Vector vz( d.Turn.zx, d.Turn.zy, d.Turn.zz );
	l.PutVector( "TX", vx );
	l.PutVector( "TY", vy );
	l.PutVector( "TZ", vz );
	l.PutVector( "S", d.Shift );
}

struct SSPred
{
	int type;
	int beg;
	int length;
	double power;
	int chain;
};

bool AddSSPred( vector<TemplateStr>& templ, vector<string>& qseq, const char *pname, double powerLimit )
{
	FILE *pfile = fopen( pname, "rt" );
	multimap<double,SSPred> pmap;
	int chainc = 0;
	do
	{
		Line l;
		if ( l.Read( pfile ) != eOK ) break;
		if ( l.GetType() == "SSPRED" )
		{
			SSPred ssp;
			ssp.type = l.GetInt( "TYPE" );
			ssp.beg = l.GetInt( "FIRST" );
			ssp.length = l.GetInt( "SIZE" );
			ssp.power = l.GetDouble( "VALUE" );
			ssp.chain = 0;
			if ( ssp.type > 0 )
			{
				pmap.insert( pair<double,SSPred>( ssp.power, ssp ) );
				if ( pmap.size() > MaxNumSS ) pmap.erase( pmap.begin() );
			}
		}
		else if ( l.GetType() == "CHAIN" )
		{
			Line ll;
			if ( l.GetFirst( ll ) ) do
			{
				SSPred ssp;
				ssp.type = ll.GetInt( "TYPE" );
				ssp.beg = ll.GetInt( "FIRST" );
				ssp.length = ll.GetInt( "SIZE" );
				ssp.power = ll.GetDouble( "VALUE" );
				ssp.chain = chainc;
				if ( ssp.type > 0 )
				{
					pmap.insert( pair<double,SSPred>( ssp.power, ssp ) );
					if ( pmap.size() > MaxNumSS ) pmap.erase( pmap.begin() );
				}
			}
			while ( l.GetNext( ll ) );
			chainc++;
		}
	}
	while ( 1 );
	fclose( pfile );
	for ( multimap<double,SSPred>::reverse_iterator it = pmap.rbegin(); it != pmap.rend(); it++ )
	{
		templ.resize( templ.size() + 1 );
		templ.back().sequence = qseq[ it->second.chain ].substr( it->second.beg, it->second.length );
		templ.back().ss.append( it->second.length, ( it->second.type == 1 ) ? 'H' : 'E' );
		templ.back().order = templ.size() + qseq.size() - 1;
		templ.back().type = 1;
		templ.back().id = -1;
		templ.back().chain = it->second.chain;
		templ.back().bb.resize( templ.back().sequence.size() );
		if ( it->second.type == 2 )
		{
			BuildBeta( templ.back().bb, it->second.length );
		}
		else
		{
			BuildAlpha( templ.back().bb, it->second.length );
		}
		templ.back().cluster.resize( 1 );
		SCluster& cluster = templ.back().cluster[ 0 ];
		cluster.h.resize( 1 );
		cluster.h[0].beg.resize( 2 );
		cluster.h[0].x = cluster.h[0].beg[0] = it->second.beg;
		cluster.h[0].y = cluster.h[0].beg[1] = 0;
		cluster.h[0].length = it->second.length;
		double power = it->second.power * 3 * sqrt( it->second.length / 6. ) ;
		cluster.sum = cluster.h[0].sum = power * sqrt( it->second.length );
		cluster.order.resize( 2 );
		cluster.order[0] = it->second.chain;
		cluster.order[1] = templ.back().order;
		cluster.power = power;
		cluster.length = it->second.length;
	}
}
		
void AddNeighbours( VBB& coord, int beg, vector<NAtom>& in, vector<NAtom>& out )
{
	double threshold = 9;
	for ( int ac = 0; ac < in.size(); ac++ )
	{
		if ( in[ac].num >= beg && in[ac].num < beg + coord.size() ) continue;
		bool find = false;
		for ( int cc = 0; cc < coord.size(); cc++ )
		{
			for ( int acc = 0; acc < 4; acc++ )
			{
				if ( ( in[ac].coord - coord[cc].v[acc] ).norm() < threshold )
				{
					find = true;
					break;
				}
			}
			if ( find ) break;
		}
		if ( find ) out.push_back( in[ac] );
	}
}
		
bool ProcessProfile( const char *sname, const char *iname, const char *pname, const char *oname )
{
	NormalizeWM( mmean, mmdisp );
	FILE *sfile = fopen( sname, "rt" );
	if ( !sfile ) return false;
	vector<StatSequence> seq;
	vector<string> qseq;
	seq.resize( 1 );
	qseq.resize( 1 );
	int chainc = 0;
	do
	{
		Line l;
		if ( l.Read( sfile ) != eOK ) break;
		if ( l.GetType() == "RES" )
		{
			StatLetter letter;
			letter.pos = l.GetInt( "POS" );
			if ( seq[0].size() && seq[0].back().pos + 1 != letter.pos ) break;
			letter.letter = l.GetString( "LETTER" )[0];
			qseq[0].append( 1, letter.letter );
			for ( int lc = 0; lc < 26; lc++ )
			{
				char curLetter = lc + 'A';
				if ( !strchr( aaLetters, curLetter ) ) continue;
				letter.distr[lc] = l.GetDouble( string( 1, curLetter ) );
			}
			seq[0].push_back( letter );
			chainc = 1;
		}
		else if ( l.GetType() == "CHAIN" )
		{
			if ( seq.size() < chainc + 1 )
			{
				seq.resize( chainc + 1 );
				qseq.resize( chainc + 1 );
			}
			Line ll;
			if ( l.GetFirst( ll ) ) do
			{
				StatLetter letter;
				letter.pos = ll.GetInt( "POS" );
				if ( seq[chainc].size() && seq[chainc].back().pos + 1 != letter.pos ) break;
				letter.letter = ll.GetString( "LETTER" )[0];
				qseq[chainc].append( 1, letter.letter );
				for ( int lc = 0; lc < 26; lc++ )
				{
					char curLetter = lc + 'A';
					if ( !strchr( aaLetters, curLetter ) ) continue;
					letter.distr[lc] = ll.GetDouble( string( 1, curLetter ) );
				}
				seq[chainc].push_back( letter );
			}
			while ( l.GetNext( ll ) );
			chainc++;
		}
	}
	while ( 1 );
	fclose( sfile );
	FILE *ifile = fopen( iname, "rt" );
	vector<TemplateStr> templ;
	int ocnt = chainc;
	SClusterSet scset;
	do
	{
		Line l;
		if ( l.Read( ifile ) != eOK ) break;
		if ( l.GetString( "SEQ" ).size() == 0 ) continue;
		TemplateStr ts;
		ts.sequence = l.GetString( "SEQ" );
		ts.ss = l.GetString( "SS" );
		ts.chain = l.GetInt( "CHAIN" );
		ts.id = l.GetInt( "ID" );
		if ( ts.chain >= seq.size() ) continue;
		if ( !ts.sequence.size() ) continue;
		templ.push_back( ts );
		Line ll;
		templ.back().bb.resize( templ.back().sequence.size() );
		templ.back().order = ocnt++;
		templ.back().type = 0;
		if ( l.GetFirst( ll ) ) do
		{
			if ( ll.GetType() == "RES" )
			{
				Backbone& bb = templ.back().bb[ ll.GetInt( "ORDER" ) ];
				bb.v[0] = ll.GetVector( "N" );
				bb.v[1] = ll.GetVector( "CA" );
				bb.v[2] = ll.GetVector( "C" );
				bb.v[3] = ll.GetVector( "O" );
			}
			else if ( ll.GetType() == "ATOM" )
			{
				NAtom na;
				strncpy( na.res, ll.GetString( "RES" ).data(), 2 );
				na.coord = ll.GetVector( "COORD" );
				strncpy( na.name, ll.GetString( "ATOM" ).data(), 4 );
				na.num = ll.GetInt( "ORDER" );
				templ.back().natom.push_back( na );
			}
		}
		while ( l.GetNext( ll ) );
		ProcessSeqDothelix( seq[ templ.back().chain ], templ.back() );
	}
	while ( 1 );
	fclose( ifile );
	cMinPower.resize( chainc );
	for ( int cc = 0; cc < chainc; cc++ )
	{
		multiset<double> pset;
		for ( int tc = 0; tc < templ.size(); tc++ )
		{
			if ( templ[tc].chain != cc ) continue;
			if( templ[tc].cluster.size() ) pset.insert( templ[tc].cluster[0].power );
			if ( pset.size() > MaxNumPClusters ) pset.erase( pset.begin() );
		}
		cMinPower[cc] = (*(pset.begin()));
		printf( "min power of cluster %g size %d\n", cMinPower[cc], ocnt );
	}
	AddSSPred( templ, qseq, pname, max( cMinPower[0], 4. ) );
	for ( int tc1 = 0; tc1 < templ.size(); tc1++ )
	{
		for ( int cc = 0; cc < templ[tc1].cluster.size(); cc++ )
		{
			if ( templ[tc1].cluster[cc].power >= cMinPower[ templ[tc1].cluster[cc].order[0] ] )
				scset.cluster.push_back( templ[tc1].cluster[cc] );
		}
		for ( int tc2 = 0; tc2 < tc1; tc2++ )
		{
			if ( templ[tc1].type == 1 && templ[tc2].type == 1 ) continue;
			if ( templ[tc1].chain != templ[tc2].chain ) continue;
			CompareTemplates( templ[tc1], templ[tc2], scset );
		}
	}
	printf( "%d initial clusters\n", scset.cluster.size() );
	scset.size = templ.size() + 1;
	SClusterSet clset;
	if ( ocnt > chainc ) 
	{
		connect_pair_clusters = false;
		FindSetClusters( scset, clset, 0 );
	}
	else 
	{
		for ( int tc = 0; tc < templ.size(); tc++ )
		{
			clset.cluster.push_back( templ[tc].cluster[0] );
		}
	}
	map<IPair,displacement> oomap;
	FILE *ofile = fopen( oname, "wt" );
	Line lh( "HEADER" );
	for ( int cc = 0; cc < qseq.size(); cc++ )
	{
		Line ll( "SEQ" );
		ll.PutInt( "CHAIN", cc );
		ll.PutString( "SEQ", qseq[cc] );
		lh.PutLine( ll );
	}
	lh.Write( ofile );
	int ccnt = 0;
	int accnt = 0;
	int fcnt = 0;
	for ( list<SCluster>::iterator it = clset.cluster.begin(); it != clset.cluster.end(); it++, ccnt++ )
	{
		vector<H3> hom3;
		vector<string> hom3ss;
		SCluster& ccl = *it;
		if ( ccl.order[0] >= chainc ) continue;
		printf( "cluster %d power %g width %d size %d\n", ccnt, ccl.power, ccl.order.size(), ccl.h.size() );
		for ( int oc = 0; oc < ccl.order.size(); oc++ ) printf( "%d ", ccl.order[oc] ); printf( "\n" );
		Line l0( "CLUSTER" );
		l0.PutInt( "ORDER", ccnt );
		l0.Write( ofile );
		int tc1, tc2;
		map<int,displacement> mmap;
		for ( tc1 = 0; tc1 < templ.size(); tc1++ )
		{
			displacement displ;
			if ( templ[tc1].order != ccl.order[1] ) continue;
			tc2 = tc1;
			int ord1 = -1;
			for ( int oc = 0; oc < ccl.order.size(); oc++ )
			{
				if ( ccl.order[oc] == templ[tc1].order ) 
				{
					ord1 = oc;
					break;
				}
			}
			for ( int hc = 0; hc < ccl.h.size(); hc++ )
			{
				Homology& h = ccl.h[hc];
				if ( h.beg[0] == -1 || h.beg[ord1] == -1 ) continue;
				H3Pos hpos( PosPair( h.beg[0], h.length ), PosPair( h.beg[ord1], h.length ) );
				H3 hom;
				hom.coord.insert( hom.coord.begin(), templ[tc1].bb.begin() + h.beg[ord1], templ[tc1].bb.begin() + h.beg[ord1] + h.length );
				hom.power = h.sum / sqrt( h.length );
				hom.pos = hpos;
				AddNeighbours( hom.coord, h.beg[ord1], templ[tc1].natom, hom.natom );
				hom3.push_back( hom );
				hom3ss.push_back( templ[tc1].ss.substr( h.beg[ord1], h.length ) );
			}
			mmap[tc1] = unary_displacement();
			int clcnt = 0;
			for ( int tc3 = tc2; tc3 < templ.size(); tc3++ )
			{
				if ( tc3 == tc1 ) continue;
				if ( !CompareTemplates2( templ[tc1], templ[tc3], ccl, displ ) ) continue;
				int ord2 = -1;
				for ( int oc = 0; oc < ccl.order.size(); oc++ )
				{
					if ( ccl.order[oc] == templ[tc3].order ) 
					{
						ord2 = oc;
						break;
					}
				}
				VBB ncoord;
				ncoord.resize( templ[tc3].bb.size() );
				vector<NAtom> natom;
				DoTransform( templ[tc3], displ, ncoord, natom );
				for ( int hc = 0; hc < ccl.h.size(); hc++ )
				{
					Homology& h = ccl.h[hc];
					if ( h.beg[0] == -1 || h.beg[ord2] == -1 ) continue;
					h.sum = 0;
					for ( int lc = 0; lc < h.length; lc++ )
					{
						char cletter = templ[tc3].sequence[ h.beg[ord2] + lc ];
						double v = 0;
						if ( cletter != 'U' && cletter != '-' ) v = seq[ ccl.order[0] ][ h.beg[0] + lc ].distr[ cletter - 'A' ];
						h.sum += v;
					}						
					H3Pos hpos( PosPair( h.beg[0], h.length ), PosPair( h.beg[ord2], h.length ) );
					H3 hom;
					hom.coord.insert( hom.coord.begin(), ncoord.begin() + h.beg[ord2], ncoord.begin() + h.beg[ord2] + h.length );
					hom.power = h.sum / sqrt( h.length );
					hom.pos = hpos;
					AddNeighbours( hom.coord, h.beg[ord2], natom, hom.natom );
					hom3.push_back( hom );
					hom3ss.push_back( templ[tc3].ss.substr( h.beg[ord2], h.length ) );
				}
				mmap[tc3] = displ;
				if ( clcnt++ > MaxNumInCluster ) break;
			}
			for ( int tc3 = 0; tc3 < templ.size(); tc3++ )
			{
				if ( templ[tc3].id == templ[tc1].id && templ[tc3].chain != templ[tc1].chain && templ[tc1].id != -1 ) 
				{
					mmap[tc3] = unary_displacement();
				}
			}
			break;
		}
		bool hasHom = false;
		for ( int hc = 0; hc < hom3.size(); hc++ )
		{
			
			H3& hom = hom3[hc];
			if ( hom.power < 0.1 ) continue;
			hasHom = true;
			Line l( "H3" );
			l.PutInt( "B", hom.pos.first.first );
			l.PutInt( "L", hom.pos.first.second );
			l.PutDouble( "P", hom.power );
			l.PutInt( "CL", ccnt );
			l.PutInt( "ORDER", fcnt++ );
			l.PutInt( "CHAIN", ccl.order[0] );
			l.PutString( "SS", hom3ss[hc] );
			for ( int cc = 0; cc < hom.coord.size(); cc++ )
			{
				Backbone& bb = hom.coord[cc];
				Line ll( "BB" );
				ll.PutInt( "N", cc );
				ll.PutVector( "A1", bb.v[0] );
				ll.PutVector( "A2", bb.v[1] );
				ll.PutVector( "A3", bb.v[2] );
				ll.PutVector( "A4", bb.v[3] );
				l.PutLine( ll );
			}
			for ( int ac = 0; ac < hom.natom.size(); ac++ )
			{
				NAtom& na = hom.natom[ac];
				Line ll( "ATOM" );
				ll.PutInt( "NUM", na.num );
				ll.PutString( "NAME", na.name );
				ll.PutString( "RES", na.res );
				ll.PutVector( "COORD", na.coord );
				l.PutLine( ll );
			}
			l.Write( ofile );
		}
		if ( hasHom )
		{
			for ( int clc = 0; clc < ccnt; clc++ )
			{
				set<int> dpoints;
				for ( map<int,displacement>::iterator it_d = mmap.begin(); it_d != mmap.end(); it_d++ )
				{
					map<IPair,displacement>::iterator it_f = oomap.find( IPair( clc, it_d->first ) );
					if ( it_f != oomap.end() )
					{
						displacement rd = multiply( it_f->second, reverse( it_d->second ) );
						if ( dpoints.find( int( rd.Turn.xx * 1000 ) ) != dpoints.end() ) continue;
						dpoints.insert( int( rd.Turn.xx * 1000 ) );
						Line ld( "DISPL" );
						ld.PutInt( "O1", clc );
						ld.PutInt( "O2", ccnt );
						ld.PutInt( "N1", -1 );
						ld.PutInt( "N2", -1 );
						WriteDisplacement( rd, ld );
						ld.Write( ofile );
					}
				}
			}
			for ( map<int,displacement>::iterator it_d = mmap.begin(); it_d != mmap.end(); it_d++ )
			{
				oomap[ IPair( ccnt, it_d->first ) ] = it_d->second;
			}
			accnt++;
			if ( accnt > MaxNumClusters ) break;
		}
	}
	fclose( ofile );
	return true;
}

int main( int argc, char **argv )
{
	if ( argc < 4 )
	{
		printf( "arguments: profile pdb sspred output_name\n" );
		return 1;
	}
	ProcessProfile( argv[1], argv[2], argv[3], argv[4] );
	return 0;
}
			
			

