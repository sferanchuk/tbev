
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <string>
#include <vector>

#include <cgicc/Cgicc.h>

using namespace std;
using namespace cgicc;

#include <pngwriter.h>

double weightMatrix[26][26] = {
{ 6, 0, -6, -3, -2, -8, -2, -7, -5, 0, -7, -6, -5, -4, 0, -2, -4, -7, 0, -1, 0, -2, -13, 0, -8, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ -6, 0, 10, -14, -14, -13, -9, -7, -6, 0, -14, -15, -13, -11, 0, -8, -14, -8, -3, -8, 0, -6, -15, 0, -4, 0},
{ -3, 0, -14, 8, 2, -15, -3, -4, -7, 0, -4, -12, -11, 2, 0, -8, -2, -10, -4, -5, 0, -8, -15, 0, -11, 0},
{ -2, 0, -14, 2, 8, -14, -4, -5, -5, 0, -4, -9, -7, -2, 0, -5, 1, -9, -4, -6, 0, -6, -17, 0, -8, 0},
{ -8, 0, -13, -15, -14, 9, -9, -6, -2, 0, -14, -3, -4, -9, 0, -10, -13, -9, -6, -9, 0, -8, -4, 0, 2, 0},
{ -2, 0, -9, -3, -4, -9, 6, -9, -11, 0, -7, -10, -8, -3, 0, -6, -7, -9, -2, -6, 0, -5, -15, 0, -14, 0},
{ -7, 0, -7, -4, -5, -6, -9, 9, -9, 0, -6, -6, -10, 0, 0, -4, 1, -2, -6, -7, 0, -6, -7, 0, -3, 0},
{ -5, 0, -6, -7, -5, -2, -11, -9, 8, 0, -6, -1, -1, -5, 0, -8, -8, -5, -7, -2, 0, 2, -14, 0, -6, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ -7, 0, -14, -4, -4, -14, -7, -6, -6, 0, 7, -8, -2, -1, 0, -6, -3, 0, -4, -3, 0, -9, -12, 0, -9, 0},
{ -6, 0, -15, -12, -9, -3, -10, -6, -1, 0, -8, 7, 1, -7, 0, -7, -5, -8, -8, -7, 0, -2, -6, 0, -7, 0},
{ -5, 0, -13, -11, -7, -4, -8, -10, -1, 0, -2, 1, 11, -9, 0, -8, -4, -4, -5, -4, 0, -1, -13, 0, -11, 0},
{ -4, 0, -11, 2, -2, -9, -3, 0, -5, 0, -1, -7, -9, 8, 0, -6, -3, -6, 0, -2, 0, -8, -8, 0, -4, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ -2, 0, -8, -8, -5, -10, -6, -4, -8, 0, -6, -7, -8, -6, 0, 8, -3, -4, -2, -4, 0, -6, -14, 0, -13, 0},
{ -4, 0, -14, -2, 1, -13, -7, 1, -8, 0, -3, -5, -4, -3, 0, -3, 8, -2, -5, -5, 0, -7, -13, 0, -12, 0},
{ -7, 0, -8, -10, -9, -9, -9, -2, -5, 0, 0, -8, -4, -6, 0, -4, -2, 8, -3, -6, 0, -8, -2, 0, -10, 0},
{ 0, 0, -3, -4, -4, -6, -2, -6, -7, 0, -4, -8, -5, 0, 0, -2, -5, -3, 6, 0, 0, -6, -5, 0, -7, 0},
{ -1, 0, -8, -5, -6, -9, -6, -7, -2, 0, -3, -7, -4, -2, 0, -4, -5, -6, 0, 7, 0, -3, -13, 0, -6, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ -2, 0, -6, -8, -6, -8, -5, -6, 2, 0, -9, -2, -1, -8, 0, -6, -7, -8, -6, -3, 0, 7, -15, 0, -7, 0},
{ -13, 0, -15, -15, -17, -4, -15, -7, -14, 0, -12, -6, -13, -8, 0, -14, -13, -2, -5, -13, 0, -15, 13, 0, -5, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{ -8, 0, -4, -11, -8, 2, -14, -3, -6, 0, -9, -7, -11, -4, 0, -13, -12, -10, -7, -6, 0, -7, -5, 0, 10, 0},
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} };

static double aa_frequences[ 26 ] = {
	0.085, 0, 0.018, 0.055, 0.066, 0.040,              // A,B,C,D,E,F
	0.075, 0.022, 0.059, 0, 0.059, 0.094,              // G,H,I,J,K,L
	0.025, 0.044, 0, 0.052, 0.040, 0.050,              // M,N,P,Q,R
	0.069, 0.060, 0, 0.072, 0.014, 0, 0.034, 0          // S,T,V,W,Y,-
	};

static double mean;
static double mdisp;

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




/*
const char *ialign =
">4pro\n"
"-LKFAMQ-RDLG----------IF--PTQLPQYLQTEKLARTQAAAIEREF-GAQ-FAGSWI\n"
">gi|561609\n"
"ELR-AIV-QDLGCGPYF-LG--TY--DKRFPGFMAPDKLACAIVNTAGRETGGEHWLAFGWN\n"
">gi|139365\n"
"ELK-AIV-RDLGCGPYF-LG--TF--DKRFPGFMAPDKLACAIVNTAGRETGGEHWLAFGWN\n"
">gi|632388\n"
"ALQ-ILYSVNLDSQRKL-INELIL---ERFGDIQENPRVLIPKNDLISRDQELSLRLQKEEE\n"
">gi|159164\n"
"FLR-DVT-KRLATDKRF----NIF--SKPVSDYLEVIKEPMDLSTVITK-IDKHNYLT----\n";
*/
enum { t0, tAll, tConsensus, tGroupConsensus, tNotConsensus, tNotGroupConsensus, tSD, tSD0, tSD1, tSD2 };

int main( int argc, char **argv )
{
	const char *ialign;
	Cgicc cgi;
	if ( cgi.getElement("align") == cgi.getElements().end() ) return 1;
	ialign = cgi.getElement( "align" )->getValue().data();
	vector<string> names;
	vector<string> seqs;
	string seq;
	int ibeg = 0;
	while ( ibeg < strlen( ialign ) )
	{
		const char *pa = strchr( ialign + ibeg, '\n' );
		int lend;
		int dr = 0;
		if ( pa == 0 ) lend = strlen( ialign );
		else 
		{
			lend = pa - ialign;
			if ( ialign[ lend - 1 ] == '\r' ) 
			{
				lend--;
				dr = 1;
			}
		}
		if ( ialign[ ibeg ] == '>' )
		{
			if ( seqs.size() < names.size() ) seqs.push_back( seq );
			seq.clear();
			names.push_back( string( ialign, ibeg + 1, lend - ibeg - 1 ) );
		}
		else
		{
			seq.append( ialign, ibeg, lend - ibeg );
		}
		ibeg = lend + 1 + dr;
	}
	if ( seqs.size() < names.size() ) seqs.push_back( seq );
	int numseq = seqs.size();
	int seqlength = 0;
	for ( int sc = 0; sc < seqs.size(); sc++ ) seqlength = max( seqlength, int( seqs[sc].size() ) );
	for ( int sc = 0; sc < seqs.size(); sc++ ) seqs[sc].append( '-',  seqlength - seqs[sc].size() );
	int namelength = 0;
	for ( int sc = 0; sc < names.size(); sc++ ) namelength = max( namelength, int( names[sc].size() ) );

	const char *fonts[] = { "FreeMono.ttf", "FreeMonoBold.ttf", "FreeSerif.ttf", "FreeSerifBold.ttf", "FreeSerifItalic.ttf", "FreeSans.ttf", "FreeSansBold.ttf", "FreeSansOblique.ttf" };
	string fontpath = "/home/sergey/work/antigen/";
	const int fontwidths[] = { 17, 17, 16, 16, 16, 16, 16, 16 };
	
	int groups[26] = { 0, 1, 2, 3, 3, 4, 0, 9, 5, 1, 6, 5, 5, 3, 1, 7, 3, 6, 8, 8, 1, 5, 4, 1, 4, 1 };
	int gap = 1;
	int seqcolor[10][3] = { { 0x8000, 0x8000, 0 }, { 0, 0, 0 }, { 0, 0xff00, 0 }, { 0x5500, 0x5500, 0xff00 }, { 0xff00, 0, 0xff00 }, { 0xff00, 0, 0 }, { 0, 0, 0xff00 }, { 0, 0xff00, 0xff00 }, { 0x5500, 0, 0xff00 }, { 0, 0xbb00, 0x7700 } };
	int cgroups[3][26] = { { 0, 1, 2, 3, 3, 4, 0, 9, 5, 1, 6, 5, 5, 3, 1, 7, 3, 6, 8, 8, 1, 5, 4, 1, 4, 1 },
	{ 0, 1, 2, 2, 3, 3, 0, 3, 2, 1, 3, 2, 3, 2, 1, 2, 3, 3, 2, 2, 1, 2, 3, 1, 3, 1 },
	{ 0, 1, 2, 3, 3, 2, 0, 3, 2, 1, 3, 2, 2, 3, 1, 3, 3, 3, 3, 3, 1, 2, 2, 1, 3, 1 } };
	
	int up_type = tConsensus;
	int down_type = tGroupConsensus;
	int shad_type = tNotGroupConsensus;
	bool use_numbers = false;
	int use_colors = tSD1;
	int fontnum = 1;
	int lfontnum = 5;
	int nfontnum = 4;
	int cgroup = 1;
	
	
	if ( cgi.getElement("uptype") != cgi.getElements().end() ) up_type = atoi( cgi.getElement( "uptype" )->getValue().data() );
	if ( cgi.getElement("downtype") != cgi.getElements().end() ) down_type = atoi( cgi.getElement( "downtype" )->getValue().data() );
	if ( cgi.getElement("shadtype") != cgi.getElements().end() ) shad_type = atoi( cgi.getElement( "shadtype" )->getValue().data() );
	if ( cgi.getElement("coltype") != cgi.getElements().end() ) use_colors = atoi( cgi.getElement( "coltype" )->getValue().data() );
	if ( cgi.getElement("numbers") != cgi.getElements().end() ) use_numbers = cgi.getElement( "numbers" )->getValue() == "1";
	if ( cgi.getElement("font") != cgi.getElements().end() ) fontnum = atoi( cgi.getElement( "font" )->getValue().data() );
	if ( cgi.getElement("lfont") != cgi.getElements().end() ) lfontnum = atoi( cgi.getElement( "lfont" )->getValue().data() );
	if ( cgi.getElement("nfont") != cgi.getElements().end() ) nfontnum = atoi( cgi.getElement( "nfont" )->getValue().data() );
	if ( cgi.getElement("group") != cgi.getElements().end() ) cgroup = atoi( cgi.getElement( "group" )->getValue().data() );
	

	NormalizeWM( mean, mdisp );
	vector<bool> consensus;
	vector<bool> groupconsensus;
	vector<double> weight;
	consensus.resize( seqlength );
	groupconsensus.resize( seqlength );
	weight.resize( seqlength );
	for ( int lc = 0; lc < seqlength; lc++ )
	{
		bool cc = true;
		bool cgc = true;
		double sum = 0;
		double norm = 0;
		for ( int sc = 0; sc < seqs.size(); sc++ )
		{
			int code1 = toupper( seqs[sc][lc] ) - 'A';
			for ( int scc = 0; scc < sc; scc++ )
			{
				if ( seqs[sc][lc] != seqs[scc][lc] ) cc = false;
				int code2 = toupper( seqs[scc][lc] ) - 'A';
				if ( code1 >= 0 && code2 >= 0 && code1 < 26 && code2 < 26 )
				{
					if ( cgroups[cgroup][code1] != cgroups[cgroup][code2] && !( cgroup == 2 && ( cgroups[2][code1] == 0 || cgroups[2][code2] == 0 ) ) ) cgc = false;
					sum += ( weightMatrix[code1][code2] - mean ) / mdisp;
					norm = norm + 1;
				}
				else cc = false;
			}
		}
		consensus[lc] = cc;
		groupconsensus[lc] = cgc;
		weight[lc] = sum / norm;
	}
	
	int maxnamelength = 20;
	int fontsize = 18;
    int fontwidth = fontwidths[fontnum];
    int lfontwidth = fontwidths[lfontnum];
    int nfontwidth = fontwidths[nfontnum];
	int numbletters = 60;
	int numblocks = seqlength / numbletters + min( seqlength % numbletters, 1 );
	int top_margin = 20;
	int left_margin = 20;
	int seqgap = 4;
	int shadoffset = 4;
	int blockgap = 20;
	int blockheight = ( numseq + min( up_type, 1 ) + min( down_type, 1 ) ) * ( fontsize + seqgap ) + blockgap;
	int pageheight = blockheight * numblocks + top_margin;
	int lengthoffset = left_margin + ( namelength + 1 ) * lfontwidth;
	int seqoffset = (use_numbers) ? lengthoffset + 4 * nfontwidth : lengthoffset;
	string font = ( fontpath + fonts[ fontnum ] );
	string lfont = ( fontpath + fonts[ lfontnum ] );
	string nfont = ( fontpath + fonts[ nfontnum ] );
	int textcolor[3] = { 0, 0, 0 };
	int pagewidth = seqoffset + fontwidth * min( numbletters, seqlength ) + left_margin;
	
	char tmpname[50];
	strcpy( tmpname, "/tmp/XXXXXX" );
	mkstemp( tmpname );
	pngwriter image( pagewidth, pageheight, 1.0, tmpname );
	
	vector<int> begs;
	begs.resize( numseq );
	for ( int sc = 0; sc < numseq; sc++ ) begs[sc] = 0;
	int abeg = 0;
	for ( int bc = 0; bc < numblocks; bc++ )
	{
		int aend = min( seqlength, abeg + numbletters );
		if ( up_type > 0 )
		{
			int top = pageheight - top_margin - bc * blockheight - fontsize / 2;
			for ( int lc = abeg; lc < aend; lc++ )
			{
				int left = seqoffset + ( lc - abeg ) * fontwidth;
				char symbol[2] = " ";
				switch ( up_type )
				{
					case tConsensus: if ( consensus[lc] ) symbol[0] = '*'; break;
					case tGroupConsensus: if ( groupconsensus[lc] ) symbol[0] = '*'; break;
					case tNotConsensus: if ( !consensus[lc] ) symbol[0] = '*'; break;
					case tNotGroupConsensus: if ( !groupconsensus[lc] ) symbol[0] = '*'; break;
					case tSD: symbol[0] = ( weight[lc] > 0 ) ? ( weight[lc] > 1 ) ? ( weight[lc] > 2 ) ? '*' : '+' : '.' : ' '; break;
				}
				image.plot_text( (char*)(font.data()), fontsize, left, top, 0, symbol, textcolor[0], textcolor[1],  textcolor[2] );
			}
		}
		for ( int sc = 0; sc < numseq; sc++ )
		{
			int top = pageheight - top_margin - bc * blockheight - ( sc + min( up_type, 1 ) ) * ( fontsize + seqgap ) - fontsize / 2;
			image.plot_text_utf8( (char*)(lfont.data()), fontsize, left_margin, top, 0, (char*)names[sc].data(), textcolor[0], textcolor[1], textcolor[2] );
			char nbuf[5];
			sprintf( nbuf, "%-3d", begs[sc] );
			if ( use_numbers ) image.plot_text( (char*)(nfont.data()), fontsize, lengthoffset, top, 0, nbuf, textcolor[0], textcolor[1], textcolor[2] );
			for ( int lc = abeg; lc < aend; lc++ )
			{
				int left = seqoffset + ( lc - abeg ) * fontwidth;
				char ch[2] = "";
				ch[0] = seqs[sc][lc];
				int *color = textcolor;
				if ( use_colors == tAll 
					|| use_colors == tConsensus && consensus[lc] 
					|| use_colors == tGroupConsensus && groupconsensus[lc] 
					|| use_colors == tNotConsensus && !consensus[lc] 
					|| use_colors == tNotGroupConsensus && !groupconsensus[lc] 
					|| use_colors == tSD0 && weight[lc] > 0
					|| use_colors == tSD1 && weight[lc] > 1
					|| use_colors == tSD2 && weight[lc] > 2
					)
				{
					int code = ch[0] - 'A';
					if ( code >= 0 && code < 26 )
					{
						int group = groups[code];
						color = &( seqcolor[group][0] );
					}
				}
				if ( shad_type == tConsensus && consensus[lc] 
					|| shad_type == tGroupConsensus && groupconsensus[lc] 
					|| shad_type == tNotConsensus && !consensus[lc] 
					|| shad_type == tNotGroupConsensus && !groupconsensus[lc] 
					|| shad_type == tSD0 && weight[lc] > 0
					|| shad_type == tSD1 && weight[lc] > 1
					|| shad_type == tSD2 && weight[lc] > 2
					)
				{
					image.filledsquare( left, top - shadoffset, left + fontwidth, top + ( fontsize + seqgap ) - shadoffset, 0.7, 0.7, 0.7 );	
				}
				image.plot_text( (char*)(font.data()), fontsize, left, top, 0, ch, color[0], color[1],  color[2] );
				
				if ( seqs[sc][lc] != ' ' && seqs[sc][lc] != '-' ) begs[sc]++;
			}
		}
		if ( down_type > 0 )
		{
			int top = pageheight - top_margin - bc * blockheight - ( numseq + min( up_type, 1 ) ) * ( fontsize + seqgap ) - fontsize / 2;
			for ( int lc = abeg; lc < aend; lc++ )
			{
				int left = seqoffset + ( lc - abeg ) * fontwidth;
				char symbol[2] = " ";
				switch ( down_type )
				{
					case tConsensus: if ( consensus[lc] ) symbol[0] = '^'; break;
					case tGroupConsensus: if ( groupconsensus[lc] ) symbol[0] = '^'; break;
					case tNotConsensus: if ( !consensus[lc] ) symbol[0] = '*'; break;
					case tNotGroupConsensus: if ( !groupconsensus[lc] ) symbol[0] = '*'; break;
					case tSD: symbol[0] = ( weight[lc] > 0 ) ? ( weight[lc] > 1 ) ? ( weight[lc] > 2 ) ? '*' : '+' : '.' : ' '; break;
				}
				image.plot_text( (char*)(font.data()), fontsize, left, top, 0, symbol, textcolor[0], textcolor[1],  textcolor[2] );
			}
		}
		abeg = aend;
	}
	image.close();
	printf( "Content-type: image/png\n\n" );
	FILE *ifile = fopen( tmpname, "rb" );
	do
	{
		int c = fgetc( ifile );
		if ( c == EOF ) break;
		fputc( c, stdout );
	}
	while ( 1 );
	fclose( ifile );
	unlink( tmpname );
}
	
		
		
	

