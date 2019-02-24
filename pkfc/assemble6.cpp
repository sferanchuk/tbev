
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <string>
#include <map>
#include <list>
#include <set>

using namespace std;

#include "vector.h"
#include "geom-f.h"
#include "format.h"
#include "cluster.h"
#include "motifs.h"
#include "motifs3d.h"
#include "dcluster.h"
#include "params.h"
#include "chain.h"

int ChainsDocking( int beg0, int beg1, PdbStruct& chain0, PdbStruct& chain1, vector<displacement>& result, bool onechain );
int ChainsQDocking( int beg0, int beg1, PdbStruct& chain0, PdbStruct& chain1, vector<NAtom>& natom0, vector<NAtom>& natom1, vector<displacement>& result, vector<double>& rpower );
int MergeChains( string seq, int beg0, int beg1, VBB& coord0, VBB& coord1, PdbStruct& chain00, PdbStruct& chain10, bool dir, double& coeff0, double& coeff1, displacement& displ, VBB& result );
void BackboneToChain( VBB& coord, int beg, string seq, PdbStruct& chain );
double CalcEnergy( PdbStruct& chain );
double CalcEnergy( VBB& coord, int beg, string seq );
bool AddSidechains( Chain& chain, Chain& schain );
bool AddHydrogens( PdbStruct& chain, PdbStruct& hchain );
bool MoveChain( PdbStruct& chain, displacement& displ );
double DRigidMove( PdbStruct& chain0, PdbStruct& chain1, displacement& result );
double RigidMove( PdbStruct& chain0, PdbStruct& chain1, displacement& result );
void ProcessTasks();
void PrintDStruct( DStruct& ds, FILE *ofile );
void TruncateHist( vector<DHist>& hist, int beg, int end );

FILE *dfile;
int printcnt = 0;

vector<DTask> theTasks;
int ServNode = 0;
extern int TimeLimit;
extern int Time0;
vector<string> theSeq;

void ReadDisplacement( displacement& d, Line& l )
{
	Vector vx = l.GetVector( "TX" );
	Vector vy = l.GetVector( "TY" );
	Vector vz = l.GetVector( "TZ" );
	d.Turn.xx = vx.x;
	d.Turn.xy = vx.y;
	d.Turn.xz = vx.z;
	d.Turn.yx = vy.x;
	d.Turn.yy = vy.y;
	d.Turn.yz = vy.z;
	d.Turn.zx = vz.x;
	d.Turn.zy = vz.y;
	d.Turn.zz = vz.z;
	d.Shift = l.GetVector( "S" );
}

bool LoadClusters( ClusterMap& cmap, vector<string>& seq, const char *iname )
{
	FILE *ifile = fopen( iname, "rt" );
	if ( !ifile )
	{
		printf( "can't open %s\n", iname );
		return false;
	}
	do
	{
		Line l;
		if ( l.Read( ifile ) != eOK ) break;
		if ( l.GetType() == "HEADER" )
		{
			Line ll;
			if ( l.GetFirst( ll ) ) do
			{
				seq.push_back( ll.GetString( "SEQ" ) );
			}
			while ( l.GetNext( ll ) );
			theSeq = seq;
			continue;
		}
		if ( l.GetType() == "DISPL" )
		{
			int o1 = l.GetInt( "O1" );
			int o2 = l.GetInt( "O2" );
			double power = l.GetDouble( "P" );
			if ( power == 0 ) power = 80;
			DDispl d;
			DDispl rd;
			ReadDisplacement( d.d, l );
			d.o1 = l.GetInt( "N1" );
			d.o2 = l.GetInt( "N2" );
			rd.o1 = d.o2;
			rd.o2 = d.o1;
			d.power = rd.power = power;
			rd.d = reverse( d.d );
			for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end(); it++ )
			{
				for ( DSMap::iterator it_s = it->second.slist.begin(); it_s != it->second.slist.end(); it_s++ )
				{
					DStruct& ds = it_s->second;
					if ( ds.ord == o1 && ( d.o1 == -1 || ds.fragments.find( d.o1 ) != ds.fragments.end() ) ) 
						ds.orientation.insert( pair<int,DDispl>( o2, d ) );
					if ( ds.ord == o2 && ( rd.o1 == -1 || ds.fragments.find( rd.o1 ) != ds.fragments.end() ) ) 
						ds.orientation.insert( pair<int,DDispl>( o1, rd ) );
				}
			}
			continue;
		}
		DStruct ds;
		int pos = l.GetInt( "B" );
		int clNum = l.GetInt( "CL" );
		int order = l.GetInt( "ORDER" );
		Line ll;
		if ( l.GetFirst( ll ) ) do
		{
			if ( ll.GetType() == "BB" )
			{
				Backbone bb;
				bb.v[0] = ll.GetVector( "A1" );
				bb.v[1] = ll.GetVector( "A2" );
				bb.v[2] = ll.GetVector( "A3" );
				bb.v[3] = ll.GetVector( "A4" );
				ds.coord.push_back( bb );
			}
			else if ( ll.GetType() == "HIST" )
			{
				DHist h;
				h.cl1 = ll.GetInt( "CL1" );
				h.cl2 = ll.GetInt( "CL2" );
				h.beg = ll.GetInt( "BEG" );
				h.end = ll.GetInt( "END" );
				int hord = ll.GetInt( "ORD" );
				if ( ds.hist.size() <= hord ) ds.hist.resize( hord + 1 );
				ds.hist[ hord ] = h;
			}
			else if ( ll.GetType() == "ATOM" )
			{
				NAtom na;
				strncpy( na.res, ll.GetString( "RES" ).data(), 2 );
				na.coord = ll.GetVector( "COORD" );
				strncpy( na.name, ll.GetString( "NAME" ).data(), 4 );
				na.num = ll.GetInt( "NUM" );
				ds.natom.push_back( na );
			}
		}
		while ( l.GetNext( ll ) );
		DDispl ud;
		ud.o1 = ud.o2 = -1;
		ud.power = 40;
		ud.d = unary_displacement();
		for ( int udc = clNum; udc <= clNum; udc++ )
		{
			ds.orientation.insert( pair<int,DDispl>( udc, ud ) );
		}
		ds.ord = clNum;
		ds.beg = pos;
		ds.sum = l.GetDouble( "P" ) * sqrt( ds.coord.size() );
		ds.chain = l.GetInt( "CHAIN" );
		double en = CalcEnergy( ds.coord, ds.beg, theSeq[ds.chain].substr( ds.beg, ds.coord.size() ) );
		ds.quality = ds.sum + en * EFactor;
		ds.power = ds.quality / sqrt( ds.coord.size() );
		ds.fragments.insert( order );
		if ( !ds.hist.size() )
		{
			ds.hist.resize( 1 );
			ds.hist[0].cl1 = ds.hist[0].cl2 = -1;
			ds.hist[0].beg = ds.beg;
			ds.hist[0].end = ds.beg + ds.coord.size();
		}
		DCluster cluster;
		cluster.slist.insert( pair<double,DStruct>( ds.quality, ds ) );
		cluster.quality = ds.power;
		cmap[ IPair( order, order ) ] = cluster;
	}
	while ( 1 );
	return true;
}

static const char *aaLetters = "ACDEFGHIKLMNPQRSTVWY";
const char *a_name[] = {"ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY",
						"HIS", "ILE", "ACE", "LYS", "LEU", "MET", "ASN",
						"no ", "PRO", "GLN", "ARG", "SER", "THR", "UNK",
						"VAL", "TRP", "no ", "TYR", "GLX"};

static const char *atomNames[4] = { "N", "CA", "C", "O" };

static void print_atom( FILE *ofile, int ac, const char *aname, const char *aatype, int rc, int lc, char chain, const char *name, Vector coord )
{
	fprintf( ofile, 
		"ATOM   %4d  %-3.3s %-3.3s %c %3d    %8.3f%8.3f%8.3f  1.00 00.00 \n",
//		"ATOM   %4d  %-3.3s %-3.3s %c %3d     %7.3f %7.3f %7.3f  1.00 00.00      %4.4s%4d\n",
		ac, aname, aatype, chain, rc, coord.x, coord.y, coord.z );
}

void SaveChain( PdbStruct& chain, const char *fname, double quality, double sum, vector<DHist>& hist )
{
	FILE *ofile = fopen( fname, "wt" );
	fprintf( ofile, "REMARK 100 MEASURE %g SUM %g\n", quality, sum );
	for ( int hc = 0; hc < hist.size(); hc++ )
	{
		fprintf( ofile, "REMARK 101 ORD %3d CL1 %4d CL2 %4d BEG %4d END %4d\n", hc, hist[hc].cl1, hist[hc].cl2, hist[hc].beg, hist[hc].end );
	}
		
	int acc = 1;
	for ( int rc = 0; rc < chain.residue.size(); rc++ )
	{
		for ( int ac = 0; ac < chain.residue[rc].atom.size(); ac++ )
		{
			Atom& a = chain.residue[rc].atom[ac];
			print_atom( ofile, acc++, a.name, chain.residue[rc].name3, chain.residue[rc].number, 0, chain.residue[rc].chain, 0, a.coord );
		}
	}
	fclose( ofile );
}

void ProcessPair( DStruct& s1, DStruct& s2, DSList& res )
{
	vector<displacement> vdispl;
	vector<displacement> vddispl;
	vector<double> dpower;
	vector<double> ddpower;
	int resOrd = -1;
	bool needrmove = false;
	if ( s1.ord == s2.ord && ( s1.beg <= s2.beg && s1.beg + s1.coord.size() > s2.beg || s2.beg < s1.beg && s2.beg + s2.coord.size() > s1.beg ) ) return;
	PdbStruct hchain0;
	PdbStruct hchain1;
	BackboneToChain( s1.coord, s1.beg, theSeq[ s1.chain ].substr( s1.beg, s1.coord.size() ), hchain0 );
	BackboneToChain( s2.coord, s2.beg, theSeq[ s2.chain ].substr( s2.beg, s2.coord.size() ), hchain1 );
	PdbStruct schain0;
	PdbStruct schain1;
	if ( !AddSidechains( hchain0, schain0 ) || !AddSidechains( hchain1, schain1 ) ) return;
	PdbStruct chain0;
	PdbStruct chain1;
	if ( !AddHydrogens( schain0, chain0 ) || !AddHydrogens( schain1, chain1 ) ) return;
	if ( s1.orientation.find( s2.ord ) != s1.orientation.end() )
	{
		for ( multimap<int,DDispl>::iterator o_it = s1.orientation.lower_bound( s2.ord ); o_it != s1.orientation.upper_bound( s2.ord ); o_it++ )
		{
			DDispl& cd = o_it->second;
			if ( cd.o2 != -1 && s2.fragments.find( cd.o2 ) == s2.fragments.end() ) continue;
			vdispl.push_back( cd.d );
			dpower.push_back( cd.power );
		}
		resOrd = s1.ord;
	}
	else if ( s2.orientation.find( s1.ord ) != s2.orientation.end() )
	{
		for ( multimap<int,DDispl>::iterator o_it = s2.orientation.lower_bound( s1.ord ); o_it != s2.orientation.upper_bound( s1.ord ); o_it++ )
		{
			DDispl& cd = o_it->second;
			if ( cd.o2 != -1 && s1.fragments.find( cd.o2 ) == s1.fragments.end() ) continue;
			vdispl.push_back( cd.d );
			dpower.push_back( cd.power );
		}
		resOrd = s2.ord;
	}
	else
	{
		ChainsQDocking( s1.beg, s2.beg, schain0, schain1, s1.natom, s2.natom, vddispl, ddpower );
		vdispl.insert( vdispl.end(), vddispl.begin(), vddispl.end() );
		dpower.insert( dpower.end(), ddpower.begin(), ddpower.end() );
		resOrd = s1.ord;
	}
	if ( vdispl.size() )
	{
		{
			int rcnt = 0;
			int dc;
			DStruct sr;
			multimap<double,displacement> sdispl;
			for ( dc = 0; dc < vdispl.size() && sdispl.size() < NumRMove; dc++ )
			{
				displacement cdispl = vdispl[dc];
				double rres = 0;
				if ( s1.ord != s2.ord && needrmove ) 
				{
					PdbStruct chain11 = chain1;
					MoveChain( chain11, cdispl );
					displacement rdispl;
					rres = RigidMove( chain0, chain11, rdispl );
					cdispl = multiply( rdispl, cdispl );
					if ( rres == 0 ) continue;
				}
				sdispl.insert( pair<double,displacement>( dpower[dc], cdispl ) );
			}
			dc = 0;
			for ( multimap<double,displacement>::reverse_iterator it = sdispl.rbegin(); it != sdispl.rend(); it++ )
			{
				VBB ncoord;
				displacement cdispl = it->second;
				int mres;
				double coeff1, coeff2;
				vector<NAtom> natom;
				if ( resOrd == s1.ord )
				{
					mres = MergeChains( theSeq[ s1.chain ], s1.beg, s2.beg, s1.coord, s2.coord, chain0, chain1, s1.sum < s2.sum, coeff1, coeff2, cdispl, ncoord );
					natom = s1.natom;
					for ( int ac = 0; ac < s2.natom.size(); ac++ )
					{
						NAtom na = s2.natom[ac];
						na.coord = cdispl.Move( na.coord );
						natom.push_back( na );
					}
				}
				else
				{
					mres = MergeChains( theSeq[ s1.chain ], s2.beg, s1.beg, s2.coord, s1.coord, chain1, chain0, s2.sum < s1.sum, coeff2, coeff1, cdispl, ncoord );
					natom = s2.natom;
					for ( int ac = 0; ac < s1.natom.size(); ac++ )
					{
						NAtom na = s1.natom[ac];
						na.coord = cdispl.Move( na.coord );
						natom.push_back( na );
					}
				}
				if ( mres != 1 ) continue;
				DStruct sr;
				sr.coord = ncoord;
				sr.natom = natom;
				sr.ord = resOrd;
				sr.beg = min( s1.beg, s2.beg );
				if ( s1.beg < s2.beg )
				{
					vector<DHist> h1 = s1.hist;
					vector<DHist> h2 = s2.hist;
					TruncateHist( h1, s1.beg, s1.beg + s1.coord.size() * ( ( resOrd == s1.ord ) ? coeff1 : coeff2 ) );
					TruncateHist( h2, s2.beg + s2.coord.size() * ( 1. - ( ( resOrd == s1.ord ) ? coeff2 : coeff1 ) ), s2.beg + s2.coord.size() );
					for ( int hc = 0; hc < h2.size(); hc++ )
					{
						if ( h2[hc].cl1 != -1 )
						{
							h2[hc].cl1 += h1.size();
							h2[hc].cl2 += h1.size();
						}
					}
					sr.hist = h1;
					sr.hist.insert( sr.hist.end(), h2.begin(), h2.end() );
					sr.hist.resize( sr.hist.size() + 1 );
					sr.hist.back().beg = sr.hist.back().end = -1;
					sr.hist.back().cl1 = h1.size() - 1;
					sr.hist.back().cl2 = sr.hist.size() - 2;
				}
				else
				{
					vector<DHist> h1 = s1.hist;
					vector<DHist> h2 = s2.hist;
					TruncateHist( h2, s2.beg, s2.beg + s2.coord.size() * ( ( resOrd == s1.ord ) ? coeff2 : coeff1 ) );
					TruncateHist( h1, s1.beg + s1.coord.size() * ( 1. - ( ( resOrd == s1.ord ) ? coeff1 : coeff2 ) ), s1.beg + s1.coord.size() );
					for ( int hc = 0; hc < h1.size(); hc++ )
					{
						if ( h1[hc].cl1 != -1 )
						{
							h1[hc].cl1 += h2.size();
							h1[hc].cl2 += h2.size();
						}
					}
					sr.hist = h2;
					sr.hist.insert( sr.hist.end(), h1.begin(), h1.end() );
					sr.hist.resize( sr.hist.size() + 1 );
					sr.hist.back().beg = sr.hist.back().end = -1;
					sr.hist.back().cl1 = h2.size() - 1;
					sr.hist.back().cl2 = sr.hist.size() - 2;
				}
				sr.chain = s1.chain;
				double en = CalcEnergy( sr.coord, sr.beg, theSeq[ sr.chain ].substr( sr.beg, sr.coord.size() ) );
				if ( en != en || en == 0 ) continue;
				if ( !needrmove ) printf( "== %d %d %d %g\n", s1.ord, s2.ord, sdispl.size(), en );
				if ( resOrd == s1.ord )
				{
					sr.orientation.insert( s1.orientation.begin(), s1.orientation.end() );
					for ( multimap<int,DDispl>::iterator o_it = s2.orientation.begin(); o_it != s2.orientation.end(); o_it++ )
					{
						if ( sr.orientation.find( o_it->first ) != sr.orientation.end() ) continue;
						DDispl nd = o_it->second;
						nd.d = multiply( cdispl, nd.d );
						sr.orientation.insert( pair<int,DDispl>( o_it->first, nd ) );
					}
				}
				else
				{
					sr.orientation.insert( s2.orientation.begin(), s2.orientation.end() );
					for ( multimap<int,DDispl>::iterator o_it = s1.orientation.begin(); o_it != s1.orientation.end(); o_it++ )
					{
						if ( sr.orientation.find( o_it->first ) != sr.orientation.end() ) continue;
						DDispl nd = o_it->second;
						nd.d = multiply( cdispl, nd.d );
						sr.orientation.insert( pair<int,DDispl>( o_it->first, nd ) );
					}
				}
				sr.sum = s1.sum * coeff1 + s2.sum * coeff2 - 0.1 * ( sr.coord.size() - coeff1 * s1.coord.size() - coeff2 * s2.coord.size() );
				//if ( needrmove ) sr.sum -= DockingPenalty * sqrt( sr.coord.size() );
				//else sr.sum -= it->first * sqrt( sr.coord.size() );
				sr.quality = sr.sum + en * EFactor + log( 1. - exp( - it->first ) ) * sqrt( s1.coord.size() + s2.coord.size() );
				sr.power = sr.quality / sqrt( sr.coord.size() );
				sr.lquality = sr.quality * sqrt( sr.coord.size() );
				sr.fragments.insert( s1.fragments.begin(), s1.fragments.end() );
				sr.fragments.insert( s2.fragments.begin(), s2.fragments.end() );
				res.push_back( sr );
				rcnt++;
				if ( rcnt > MaxDocking ) break;
			}
		}
	}
}

void DockingPair( DStruct& s1, DStruct& s2, multimap<double,PdbStruct>& res )
{
	vector<displacement> vdispl;
	vector<displacement> vddispl;
	vector<double> dpower;
	vector<double> ddpower;
	int resOrd = -1;
	bool needrmove = false;
	if ( s1.ord == s2.ord && ( s1.beg <= s2.beg && s1.beg + s1.coord.size() > s2.beg || s2.beg < s1.beg && s2.beg + s2.coord.size() > s1.beg ) ) return;
	PdbStruct hchain0;
	PdbStruct hchain1;
	printf( "sd chains %d %d\n", s1.chain, s2.chain );
	BackboneToChain( s1.coord, s1.beg, theSeq[ s1.chain ].substr( s1.beg, s1.coord.size() ), hchain0 );
	BackboneToChain( s2.coord, s2.beg, theSeq[ s2.chain ].substr( s2.beg, s2.coord.size() ), hchain1 );
	PdbStruct schain0;
	PdbStruct schain1;
	if ( !AddSidechains( hchain0, schain0 ) || !AddSidechains( hchain1, schain1 ) ) return;
	PdbStruct chain0;
	PdbStruct chain1;
	if ( !AddHydrogens( schain0, chain0 ) || !AddHydrogens( schain1, chain1 ) ) return;
	for ( int rc = 0; rc < chain0.residue.size(); rc++ ) chain0.residue[rc].chain = s1.chain + 'A';
	for ( int rc = 0; rc < chain1.residue.size(); rc++ ) 
	{
		chain1.residue[rc].chain = s2.chain + 'A';
		chain1.residue[rc].number += theSeq[ s1.chain ].size();
	}
	if ( s1.orientation.find( s2.ord ) != s1.orientation.end() )
	{
		for ( multimap<int,DDispl>::iterator o_it = s1.orientation.lower_bound( s2.ord ); o_it != s1.orientation.upper_bound( s2.ord ); o_it++ )
		{
			DDispl& cd = o_it->second;
			if ( cd.o2 != -1 && s2.fragments.find( cd.o2 ) == s2.fragments.end() ) continue;
			vdispl.push_back( cd.d );
			dpower.push_back( cd.power );
		}
		resOrd = s1.ord;
	}
	else if ( s2.orientation.find( s1.ord ) != s2.orientation.end() )
	{
		for ( multimap<int,DDispl>::iterator o_it = s2.orientation.lower_bound( s1.ord ); o_it != s2.orientation.upper_bound( s1.ord ); o_it++ )
		{
			DDispl& cd = o_it->second;
			if ( cd.o2 != -1 && s1.fragments.find( cd.o2 ) == s1.fragments.end() ) continue;
			vdispl.push_back( cd.d );
			dpower.push_back( cd.power );
		}
		resOrd = s2.ord;
	}
	else
	{
		ChainsQDocking( s1.beg, s2.beg, schain0, schain1, s1.natom, s2.natom, vddispl, ddpower );
		vdispl.insert( vdispl.end(), vddispl.begin(), vddispl.end() );
		dpower.insert( dpower.end(), ddpower.begin(), ddpower.end() );
		resOrd = s1.ord;
	}
	printf( "cd vdispl %d\n", vdispl.size() );
	if ( vdispl.size() )
	{
		{
			int rcnt = 0;
			int dc;
			DStruct sr;
			multimap<double,displacement> sdispl;
			for ( dc = 0; dc < vdispl.size() && sdispl.size() < NumRMove; dc++ )
			{
				displacement cdispl = vdispl[dc];
				double rres = 0;
				if ( s1.ord != s2.ord && needrmove ) 
				{
					PdbStruct chain11 = chain1;
					MoveChain( chain11, cdispl );
					displacement rdispl;
					rres = RigidMove( chain0, chain11, rdispl );
					cdispl = multiply( rdispl, cdispl );
					if ( rres == 0 ) continue;
				}
				sdispl.insert( pair<double,displacement>( dpower[dc], cdispl ) );
			}
			dc = 0;
			printf( "cd sdispl %d\n", sdispl.size() );
			for ( multimap<double,displacement>::reverse_iterator it = sdispl.rbegin(); it != sdispl.rend(); it++ )
			{
				VBB ncoord;
				displacement cdispl = it->second;
				PdbStruct chain11 = chain1;
				MoveChain( chain11, cdispl );
				chain11.residue.insert( chain11.residue.begin(), chain0.residue.begin(), chain0.residue.end() );
				double en = CalcEnergy( chain11 );
				double quality = s1.quality + s2.quality + log( 1. - exp( - it->first ) ) * sqrt( s1.coord.size() + s2.coord.size() ) + en * EFactor;
				res.insert( pair<double,PdbStruct>( quality, chain11 ) );
				rcnt++;
				if ( rcnt > MaxDocking ) break;
			}
		}
	}
}

/*
void CompareClusters( string& seq, DCluster& c1, DCluster& c2, DCluster& res )
{
	res.quality = 0;
	for ( DSList::iterator it1 = c1.slist.begin(); it1 != c1.slist.end(); it1++ )
	{
		DStruct& s1 = *it1;
		for ( DSList::iterator it2 = c2.slist.begin(); it2 != c2.slist.end(); it2++ )
		{
			DStruct& s2 = *it2;
			vector<displacement> vdispl;
			int resOrd = -1;
			if ( s1.orientation.find( s2.ord ) != s1.orientation.end() )
			{
				for ( multimap<int,displacement>::iterator o_it = s1.orientation.lower_bound( s2.ord ); o_it != s1.orientation.upper_bound( s2.ord ); o_it++ )
				{
					vdispl.push_back( o_it->second );
				}
				resOrd = s1.ord;
			}
			else if ( s2.orientation.find( s1.ord ) != s2.orientation.end() )
			{
				for ( multimap<int,displacement>::iterator o_it = s2.orientation.lower_bound( s1.ord ); o_it != s2.orientation.upper_bound( s1.ord ); o_it++ )
				{
					vdispl.push_back( o_it->second );
				}
				resOrd = s2.ord;
			}
			else
			{
				ChainsDocking( seq, s1.beg, s2.beg, s1.coord, s2.coord, vdispl );
				resOrd = s1.ord;
			}
			int rcnt = 0;
			int dc;
			DStruct sr;
			for ( dc = 0; dc < vdispl.size(); )
			{
				VBB ncoord;
				displacement cdispl = vdispl[dc];
				double mres;
				if ( resOrd == s1.ord )
				{
					mres = MergeChains( seq, s1.beg, s2.beg, s1.coord, s2.coord, cdispl, ncoord );
				}
				else
				{
					mres = MergeChains( seq, s2.beg, s1.beg, s2.coord, s1.coord, cdispl, ncoord );
				}
				if ( mres == 0 )
				{
					vdispl.erase( vdispl.begin() + dc );
				}
				else
				{
					fprintf( dfile, "++++ %d %d size %d %g %g\n", s1.ord, s2.ord, vdispl.size(), s1.quality + s2.quality, mres );
					vdispl[dc] = cdispl;
					dc++;
					DStruct sr;
					sr.coord = ncoord;
					sr.ord = resOrd;
					sr.beg = min( s1.beg, s2.beg );
					if ( resOrd == s1.ord )
					{
						sr.orientation.insert( s1.orientation.begin(), s1.orientation.end() );
						for ( multimap<int,displacement>::iterator o_it = s2.orientation.begin(); o_it != s2.orientation.end(); o_it++ )
						{
							if ( sr.orientation.find( o_it->first ) != sr.orientation.end() ) continue;
							sr.orientation.insert( pair<int,displacement>( o_it->first, multiply( cdispl, o_it->second ) ) );
						}
					}
					else
					{
						sr.orientation.insert( s2.orientation.begin(), s2.orientation.end() );
						for ( multimap<int,displacement>::iterator o_it = s1.orientation.begin(); o_it != s1.orientation.end(); o_it++ )
						{
							if ( sr.orientation.find( o_it->first ) != sr.orientation.end() ) continue;
							sr.orientation.insert( pair<int,displacement>( o_it->first, multiply( cdispl, o_it->second ) ) );
						}
					}
					sr.quality = s1.quality + s2.quality - mres;
					if ( res.quality < sr.quality ) res.quality = sr.quality;
					res.slist.push_back( sr );
					rcnt++;
					if ( rcnt > MaxDocking ) break;
				}
			}
		}
	}
}
*/

void ClusterizeParallel( ClusterMap& cmap )
{
	int nclusters = cmap.size();
	for ( int cc = 0; cc < nclusters; cc++ )
	{
		DCluster& c1 = cmap[ IPair( cc, cc ) ];
		for ( int cc1 = 0; cc1 < cc; cc1++ )
		{
			DCluster& c2 = cmap[ IPair( cc1, cc1 ) ];
			for ( DSMap::iterator it1 = c1.slist.begin(); it1 != c1.slist.end(); it1++ )
			{
				DStruct& s1 = it1->second;
				for ( DSMap::iterator it2 = c2.slist.begin(); it2 != c2.slist.end(); it2++ )
				{
					DStruct& s2 = it2->second;
					DTask task;
					task.input1 = s1;
					task.input2 = s2;
					task.b1 = cc;
					task.b2 = cc1;
					task.posted = false;
					task.completed = false;
					theTasks.push_back( task );
				}
			}
		}
	}
	ProcessTasks();
	for ( int tc = 0; tc < theTasks.size(); tc++ )
	{
		if ( !theTasks[tc].output.size() ) continue;
		IPair pos( theTasks[tc].b1, theTasks[tc].b2 );
		if ( cmap.find( pos ) == cmap.end() ) cmap[ pos ].quality = 0;
		for ( DSList::iterator it = theTasks[tc].output.begin(); it != theTasks[tc].output.end(); it++ )
		{
			cmap[pos].slist.insert( pair<double,DStruct>( it->quality, *it ) );
			if ( cmap[pos].slist.size() > CashSize ) cmap[pos].slist.erase( cmap[pos].slist.begin() );
			cmap[pos].quality = fmax( cmap[pos].quality, it->power );
		}
	}
	for ( int ic = 0; ic < nclusters; ic++ )
	{
		theTasks.clear();
		int b1;
		int b2;
		double bestquality = 0;
		printf( "ic %d\n", ic );
		for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end(); it++ )
		{
			if ( it->first.first == it->first.second ) continue;
			if ( it->second.quality > bestquality )
			{
				bestquality = it->second.quality;
				b1 = it->first.first;
				b2 = it->first.second;
			}
		}
		if ( bestquality == 0 ) break;
		printf( "clusterizing %d %d at %g\n", b1, b2, bestquality );
		DCluster pres = cmap[ IPair( b1, b2 ) ];
		for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end();  )
		{
			if ( it->first.first == b2 || it->first.second == b2 )
			{
				cmap.erase( it++ );
				continue;
			}
			else if ( it->first.first == b1 && it->first.second == b1 )
			{
				IPair ip = it->first;
				it++;
				cmap[ ip ] = pres;
			}
			else if ( it->first.first == it->first.second )
			{
				for ( DSMap::iterator it1 = it->second.slist.begin(); it1 != it->second.slist.end(); it1++ )
				{
					DStruct& s1 = it1->second;
					for ( DSMap::iterator it2 = pres.slist.begin(); it2 != pres.slist.end(); it2++ )
					{
						DStruct& s2 = it2->second;
						DTask task;
						task.input1 = s1;
						task.input2 = s2;
						task.b1 = max( b1, it->first.second );
						task.b2 = min( b1, it->first.second );
						task.posted = false;
						task.completed = false;
						theTasks.push_back( task );
					}
				}
				it++;
			}
			else it++;
		}
		ProcessTasks();
		for ( int tc = 0; tc < theTasks.size(); tc++ )
		{
			IPair pos( theTasks[tc].b1, theTasks[tc].b2 );
			cmap[pos].quality = 0;
			for ( DSList::iterator it = theTasks[tc].output.begin(); it != theTasks[tc].output.end(); it++ )
			{
				cmap[pos].slist.insert( pair<double,DStruct>( it->quality, *it ) );
				if ( cmap[pos].slist.size() > CashSize ) cmap[pos].slist.erase( cmap[pos].slist.begin() );
				cmap[pos].quality = fmax( cmap[pos].quality, it->power );
			}
		}
	}
}

/*
void Clusterize( ClusterMap& cmap, string& seq )
{
	int nclusters = cmap.size();
	for ( int cc = 0; cc < nclusters; cc++ )
	{
		DCluster& c1 = cmap[ IPair( cc, cc ) ];
		for ( int cc1 = 0; cc1 < cc; cc1++ )
		{
			DCluster& c2 = cmap[ IPair( cc1, cc1 ) ];
			DCluster res;
			CompareClusters( seq, c1, c2, res );
			if ( !res.slist.size() ) continue;
			cmap[ IPair( cc1, cc ) ] = res;
		}
	}
	for ( int ic = 0; ic < nclusters; ic++ )
	{
		int b1;
		int b2;
		double bestquality = 0;
		for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end(); it++ )
		{
			if ( it->first.first == it->first.second ) continue;
			if ( it->second.quality > bestquality )
			{
				bestquality = it->second.quality;
				b1 = it->first.first;
				b2 = it->first.second;
			}
		}
		if ( bestquality == 0 ) break;
		DCluster pres = cmap[ IPair( b1, b2 ) ];
		for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end();  )
		{
			if ( it->first.first == b2 || it->first.second == b2 )
			{
				cmap.erase( it++ );
				continue;
			}
			else if ( it->first.first == b1 && it->first.second == b1 )
			{
				IPair ip = it->first;
				it++;
				cmap[ ip ] = pres;
			}
			else if ( it->first.first == b1 || it->first.second == b1 )
			{
				DCluster nres;
				CompareClusters( seq, it->second, pres, nres );
				IPair ip = it->first;
				it++;
				if ( nres.slist.size() ) cmap[ ip ] = nres;
			}
			else it++;
		}
	}
}
*/

bool Assemble( const char *input, const char *output )
{
	ClusterMap cmap;
	vector<string> seq;
	if ( !LoadClusters( cmap, seq, input ) ) return false;
	char dname[80];
	sprintf( dname, "%s.log", output );
	dfile = fopen( dname, "wt" );
	Line lseq( "HEADER" );
	for ( int cc = 0; cc < seq.size(); cc++ )
	{
		Line ll( "SEQ" );
		ll.PutInt( "CHAIN", cc );
		ll.PutString( "SEQ", seq[cc] );
		lseq.PutLine( ll );
	}
	lseq.Write( dfile );
	for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end(); it++ )
	{
		for ( DSMap::iterator it1 = it->second.slist.begin(); it1 != it->second.slist.end(); it1++ )
		{
			PrintDStruct( it1->second, dfile );
		}
	}
	DSMap result;
	DSList rlist[2];
	for ( int cc = 0; cc < seq.size(); cc++ )
	{
		ClusterMap ccmap;
		for ( ClusterMap::iterator it = cmap.begin(); it != cmap.end(); it++ )
		{
			if ( it->second.slist.begin()->second.chain == cc ) ccmap.insert( *it );
		}
		ClusterizeParallel( ccmap );
		DSMap cresult;
		for ( ClusterMap::iterator it = ccmap.begin(); it != ccmap.end(); it++ )
		{
			for ( DSMap::iterator it1 = it->second.slist.begin(); it1 != it->second.slist.end(); it1++ )
			{
				cresult.insert( pair<double,DStruct>( it1->second.lquality, it1->second ) );
				if ( cresult.size() > 5 ) cresult.erase( cresult.begin() );
			}
		}
		result.insert( cresult.begin(), cresult.end() );
		for ( DSMap::reverse_iterator it = cresult.rbegin(); it != cresult.rend(); it++ )
		{
			DSList& clist = ( cc == 0 ) ? rlist[0] : rlist[1];
			clist.push_back( it->second );
		}
	}
	if ( seq.size() == 1 )
	{
		int ocnt = 1;
		for ( DSMap::reverse_iterator it = result.rbegin(); it != result.rend(); it++, ocnt++ )
		{
			DStruct& res = it->second;
			if ( ocnt <= 5 )
			{
				PdbStruct chain;
				BackboneToChain( res.coord, res.beg, seq[0].substr( res.beg, res.coord.size() ), chain );
				char ofmt[80];
				strcpy( ofmt, output );
				strcat( ofmt, "-m%d.pdb" );
				char obuf[100];
				sprintf( obuf, ofmt, ocnt );
				SaveChain( chain, obuf, res.quality, res.sum, res.hist );
			}
		}
	}
	else
	{
		multimap<double,PdbStruct> rstruct;
		for ( DSIt it1 = rlist[0].begin(); it1 != rlist[0].end(); it1++ )
		{
			for ( DSIt it2 = rlist[1].begin(); it2 != rlist[1].end(); it2++ )
			{
				DockingPair( *it1, *it2, rstruct );
			}
		}
		int ocnt = 1;
		vector<DHist> emptyHist;
		for ( multimap<double,PdbStruct>::reverse_iterator it = rstruct.rbegin(); it != rstruct.rend() && ocnt <= 5; it++, ocnt++ )
		{
			char ofmt[80];
			strcpy( ofmt, output );
			strcat( ofmt, "-m%d.pdb" );
			char obuf[100];
			sprintf( obuf, ofmt, ocnt );
			SaveChain( it->second, obuf, it->first, 0, emptyHist );
		}
	}
	fclose( dfile );
	return true;
}


#ifdef __MAIN__

void ProcessTasks()
{
	for ( int tc = 0; tc < theTasks.size(); tc++ )
	{
		DTask& task = theTasks[tc];
		ProcessPair( task.input1, task.input2, task.output );
	}
}

bool DeleteHist( vector<DHist>& hist, int pos )
{
	int jpos = -1;
	int kpos;
	for ( int hc = 0; hc < hist.size(); hc++ )
	{
		if ( hist[hc].cl1 == pos )
		{
			jpos = hist[hc].cl2;
			kpos = hc;
		}
		if ( hist[hc].cl2 == pos )
		{
			jpos = hist[hc].cl1;
			kpos = hc;
		}
	}
	if ( jpos == -1 ) return false;
	if ( jpos != -1 )
	{
		for ( int hc = 0; hc < hist.size(); hc++ )
		{
			if ( hist[hc].cl1 == kpos ) hist[hc].cl1 = jpos;
			if ( hist[hc].cl2 == kpos ) hist[hc].cl2 = jpos;
		}
	}
	for ( int hc = 0; hc < hist.size(); hc++ )
	{
		if ( hist[hc].cl1 > pos ) hist[hc].cl1--;
		if ( hist[hc].cl2 > pos ) hist[hc].cl2--;
	}
	hist.erase( hist.begin() + pos );
	if ( jpos != -1 )
	{
		if ( kpos > pos ) kpos--;
		for ( int hc = 0; hc < hist.size(); hc++ )
		{
			if ( hist[hc].cl1 > kpos ) hist[hc].cl1--;
			if ( hist[hc].cl2 > kpos ) hist[hc].cl2--;
		}
		hist.erase( hist.begin() + kpos );
	}
	return true;
}
	
void TruncateHist( vector<DHist>& hist, int beg, int end )
{
	do
	{
		bool find = false;
		for ( int hc = 0; hc < hist.size(); hc++ )
		{
			if ( hist[hc].cl1 == -1 && ( hist[hc].beg >= end || hist[hc].end <= beg ) )
			{
				if ( !DeleteHist( hist, hc ) )
				{
					printf( "hist error %d %d %d %d\n", hist.size(), hc, beg, end );
					exit( 1 );
				}
				find = true;
				break;
			}
		}
		if ( !find ) break;
	}
	while ( 1 );
	for ( int hc = 0; hc < hist.size(); hc++ )
	{
		if ( hist[hc].cl1 == -1 && hist[hc].beg < beg ) hist[hc].beg = beg;
		if ( hist[hc].cl1 == -1 && hist[hc].end > end ) hist[hc].end = end;
	}
}

extern int printcnt;

void PrintDStruct( DStruct& ds, FILE *ofile )
{
	Line lo( "H3" );
	lo.PutInt( "B", ds.beg );
	lo.PutInt( "L", ds.coord.size() );
	lo.PutDouble( "S", ds.sum );
	lo.PutDouble( "Q", ds.quality );
	lo.PutDouble( "P", ds.power );
	lo.PutInt( "CL", ds.ord );
	lo.PutInt( "CHAIN", ds.chain );
	lo.PutInt( "ORDER", printcnt++ );
	for ( int cc = 0; cc < ds.coord.size(); cc++ )
	{
		Backbone& bb = ds.coord[cc];
		Line ll( "BB" );
		ll.PutInt( "N", cc );
		ll.PutVector( "A1", bb.v[0] );
		ll.PutVector( "A2", bb.v[1] );
		ll.PutVector( "A3", bb.v[2] );
		ll.PutVector( "A4", bb.v[3] );
		lo.PutLine( ll );
	}
	for ( int hc = 0; hc < ds.hist.size(); hc++ )
	{
		DHist& h = ds.hist[hc];
		Line ll( "HIST" );
		ll.PutInt( "ORD", hc );
		ll.PutInt( "CL1", h.cl1 );
		ll.PutInt( "CL2", h.cl2 );
		ll.PutInt( "BEG", h.beg );
		ll.PutInt( "END", h.end );
		lo.PutLine( ll );
	}
	lo.Write( ofile );
}


int main( int argc, char **argv )
{
	if ( argc < 3 )
	{
		printf( "arguments: input_motifs output_name [time_limit]\n" );
		return 1;
	}
	if ( argc > 3 ) TimeLimit = atoi( argv[3] );
	Time0 = time( 0 );
	Assemble( argv[1], argv[2] );
	return 0;
}

#endif

		
				
			
			