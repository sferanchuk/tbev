
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
using namespace std;

#include "vector.h"
#include "geom-f.h"

struct Atom
{
	Vector coord;
	int number;
	char name[5];
};

enum { aN, aCA, aC, aO };

struct Residue
{
	vector<Atom> atom;
	Vector center;
	int number;
	char icode;
	char name3[5];
	char name1;
	int order;
	int sstype;
	int rotamer;
	double freq;
};

struct Chain
{
	vector<Residue> residue;
};

struct PDB
{ 
	vector<Chain> chain;
};

const int size1 = 5;
const int size2 = 2;
const double threshold = 0.2;

const int pdbLineSize = 100;

static const char *a_name[] = {"ALA", "ASX", "CYS", "ASP", "GLU", "PHE", "GLY",
						"HIS", "ILE", "ACE", "LYS", "LEU", "MET", "ASN",
						"no ", "PRO", "GLN", "ARG", "SER", "THR", "UNK",
						"VAL", "TRP", "no ", "TYR", "GLX"};

struct pdbLine
{
	char s[ pdbLineSize ];

	char *lineType() { return s; }
	char *atomNumber() { return s + 7; }
	char *atomType() { return s + 13; }
	char *chainName() { return s + 21; }
	char *aaNumber() { return s + 23; }
	char *aaName() { return s + 17; }
	char *x() { return s + 31; }
	char *y() { return s + 39; }
	char *z() { return s + 47; }
};

void LoadChain( PDB& pdb, const char *fname )
{
	FILE *ifile = fopen( fname, "rt" );
	pdbLine line;
	pdb.chain.resize( 1 );
	Chain *ch = &( pdb.chain.back() );
	bool first = true;
	while( ifile && fgets( line.s, pdbLineSize, ifile ) )
	{
		if ( strncmp( line.lineType(), "MODEL", 5 ) == 0 )
		{
			if ( !first )
			{
				pdb.chain.resize( pdb.chain.size() + 1 );
				ch = &( pdb.chain.back() );
			}
			else first = false;
		}
		if ( strncmp( line.lineType(), "ATOM", 4 ) ) continue;
		if ( !ch->residue.size() || ch->residue.back().number != atoi( line.aaNumber() ) )
		{
			ch->residue.resize( ch->residue.size() + 1 );
			strncpy( ch->residue.back().name3, line.aaName(), 3 );
			ch->residue.back().name3[3] = 0;
			for ( int lc = 0; lc < 26; lc++ ) if ( strncmp( line.aaName(), a_name[lc], 3 ) == 0 ) 
				ch->residue.back().name1 = lc + 'A';
			ch->residue.back().number = atoi( line.aaNumber() );
		}
		Atom a;
		memset( a.name, 0, sizeof( a.name ) );
		strncpy( a.name, line.atomType(), strcspn( line.atomType(), " " ) );
		if ( strcmp( a.name, "OXT" ) == 0 ) strcpy( a.name, "O" );
		a.coord = Vector( atof( line.x() ), atof( line.y() ), atof( line.z() ) );
		ch->residue.back().atom.push_back( a );
	}
	if ( ifile ) fclose( ifile );
}

displacement *CompareStructures( Chain& ch1, Chain& ch2 )
{
	vector<Vector> v1;
	vector<Vector> v2;
	v1.resize( ch1.residue.size() );
	v2.resize( ch2.residue.size() );
	for ( int rc = 0; rc < min( v1.size(), v2.size() ); rc++ )
	{
		v1[rc] = ch1.residue[rc].atom[aCA].coord;
		v2[rc] = ch2.residue[rc].atom[aCA].coord;
	}
	return FindDisplacement( min( v1.size(), v2.size() ), &( v1[0] ), &( v2[0] ) );
}

double AddRmsd( Chain& ch1, Chain& ch2, displacement *displ, displacement *d0, double *sum, Vector *sumv )
{
	double rv = 0;
	int cnt1 = 0;
	for ( int rc = 0; rc < min( ch1.residue.size(), ch2.residue.size() ); rc++ )
	{
		Vector v1 = d0->Move( ch1.residue[rc].atom[aCA].coord );
		Vector v2 = d0->Move( displ->Move( ch2.residue[rc].atom[aCA].coord ) );
		double d = ( v1 - v2 ).norm();
		sum[rc] += d * d;
		sumv[rc] = sumv[rc] + ( v1 - v2 );
		if ( ( rc >= 47 && rc <= 54 ) || ( rc >= 204 && rc <= 208 ) )
		{
			rv += d * d;
			cnt1++;
		}
	}
	return sqrt( rv ) / cnt1;
}

static double mean[1000] = { 0 };

void ProcessPair( PDB& pdb1, PDB& pdb2, const char *fname, const char *vfname )
{
	int size = min( pdb1.chain[0].residue.size(), pdb2.chain[0].residue.size() );
	double *sum = new double[ size ];
	Vector *sumv = new Vector[ size ];
	memset( sum, 0, size * sizeof( double ) );
	memset( sumv, 0, size * sizeof( Vector ) );
	int cnt = 0;
	int n1;
	int n2;
	double maxrmsd = 0;
	for ( int c1 = 0; c1 < pdb1.chain.size(); c1++ )
	{
		Vector v11 = pdb1.chain[c1].residue[17].atom[0].coord;
		Vector v12 = pdb1.chain[c1].residue[149].atom[0].coord;
		if ( ( v11 - v12 ).norm() > 16 ) continue;
		displacement *d0 = CompareStructures( pdb1.chain[0], pdb1.chain[c1] );
		
		for ( int c2 = 0; c2 < pdb2.chain.size(); c2++ )
		{
			Vector v21 = pdb2.chain[c2].residue[17].atom[0].coord;
			Vector v22 = pdb2.chain[c2].residue[149].atom[0].coord;
			if ( ( v21 - v22 ).norm() > 16 ) continue;

			displacement *displ = CompareStructures( pdb1.chain[c1], pdb2.chain[c2] );
			double rmsd = AddRmsd( pdb1.chain[c1], pdb2.chain[c2], displ, d0, sum, sumv );
			delete displ;
			if ( rmsd > maxrmsd )
			{
				n1 = c1;
				n2 = c2;
				maxrmsd = rmsd;
			}
			cnt++;
		}
		delete d0;
	}
	FILE *ofile = fopen( fname, "wt" );
	for ( int oc = 0; oc < size; oc++ )
	{
		fprintf( ofile, "%d %g\n", oc, sqrt( sum[oc] / cnt ) );
		if ( !vfname ) mean[oc] += sqrt( sum[oc] / cnt );
	}
	fclose( ofile );
	if ( vfname )
	{
		FILE *vfile = fopen( vfname, "wt" );
		for ( int oc = 0; oc < size; oc++ )
		{
			fprintf( ofile, "%d %g\n", oc, sumv[oc].norm() / cnt );
		}
		fclose( ofile );
	}
	delete sum;
	delete sumv;
	printf( "max rmsd: n1 = %d, n2 = %d rmsd %g\n", n1, n2, maxrmsd );
}

int main( int argc, char **argv )
{
	if ( argc < 3 ) 
	{
		printf( "usage: pdb1 pdb2\n" );
		return 1;
	}
	char fname1[50], fname2[50], fname3[50], fname4[50];
	strcpy( fname1, argv[1] );
	strtok( fname1, "." );
	strcpy( fname2, argv[2] );
	strtok( fname2, "." );
	sprintf( fname3, "%s-%s-cmp.dat", fname1, fname2 );
	sprintf( fname4, "%s-%s-vcmp.dat", fname1, fname2 );
	strcat( fname1, "-cmp.dat" );
	strcat( fname2, "-cmp.dat" );
	PDB pdb1;
	PDB pdb2;
	LoadChain( pdb1, argv[1] );
	LoadChain( pdb2, argv[2] );
	ProcessPair( pdb1, pdb1, fname1, 0 );
	ProcessPair( pdb2, pdb2, fname2, 0 );
	ProcessPair( pdb1, pdb2, fname3, fname4 );
	return 0;
}
