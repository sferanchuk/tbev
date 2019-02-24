
#include <string>
#include <vector>

using namespace std;

#include <pngwriter.h>

void ag_freq( const char *seq, int *result, double strength, double increment, int length_threshold );
int ind( char l );
extern double aa_weight[];

class gbPlotPrintout
{
 public:
  char *id;
  vector<string> letter;
  vector< vector<double> > value;
  vector< vector<double> > color;
  int beg;
  int end;
  int offset;
  double threshold;
  double maxPlot;

  gbPlotPrintout() 
  {
		beg = 0;
		end = 10000;
		offset = 0;
  }
  void DrawPage( pngwriter& image );
};

void gbPlotPrintout::DrawPage( pngwriter& image )
{
	char angstroms[10] = "5 \xc3\x85";
	int length = min( letter[0].size(), value[0].size() );
	/*
    double maxPlot = 0, minPlot = 1;
    for ( int lc = beg; lc < end; lc++ )
    {
		for ( int vc = 0; vc < value.size(); vc++ )
		{
			if ( value[vc][lc] > maxPlot ) maxPlot = value[vc][lc];
			if ( value[vc][lc] < minPlot ) minPlot = value[vc][lc];
		}
    } 
    fprintf( stderr, "%g %g %d\n", minPlot, maxPlot, int( value.size() ) );
    */
    double minPlot = 0;
	//maxPlot = printout.maxPlot;

    int w = image.getwidth();
	int h = image.getheight();
    
    // We know the graphic is 200x200. If we didn't know this,
    // we'd need to calculate it.
    float maxX = 1000;
    float maxY = 1000;
    
    // Let's have at least 50 device units margin
    float marginX = 100;
    float marginY = 100;
    
    // Get the size of the DC in pixels
    
    // Calculate a suitable scaling factor
    float scaleX=(float)(w/( maxX + 2 * marginX  ) );
    float scaleY=(float)(h/( maxY + 2 * marginY ) );
    
    // Use x or y scaling factor, whichever fits on the DC
    float actualScale = scaleX < scaleY ? scaleX : scaleY;
    
    // Calculate the position on the DC for centring the graphic
    float fPosX = (float)((w - (maxX*scaleX))/2.0);
    float fPosY = (float)((h - (maxY*scaleY))/2.0);

	int numparts = 10;
	int numletters = 10 * numparts;
	int fontsize = 20;
    int posX = 0, posY = h;
    int width = 16;
    int interval = ( h - 1  ) / ( ( end - beg ) / numletters  + 1 );
	if ( interval > 280 ) interval = 280;
    int plotHeight = interval - 10 - fontsize * letter.size();

    for ( int lc = beg; lc <= end; lc++ )
    {
    	if ( lc > beg && ( ( lc - beg ) / 10 ) % numparts != 0 && ( lc - beg )  % 10 == 0 )
    	{
    		image.line( posX, posY, posX, posY - interval, 0.6, 0.6, 0.6 );
    	}
    	if ( lc == end ) break;
    	if ( posX / width >= numletters || lc == end - 1 )
    	{
			for ( double scale = 0; scale < maxPlot; scale += 2 )
			{
				image.line( 0, posY - plotHeight + ( scale - minPlot ) / ( maxPlot - minPlot ) * plotHeight, 
					w, posY - plotHeight + ( scale - minPlot ) / ( maxPlot - minPlot ) * plotHeight, 0.85, 0.85, 0.85 );
				//image.plot_text_utf8( "./FreeMono.ttf", fontsize, 10, posY - plotHeight + ( 5 - minPlot ) / ( maxPlot - minPlot ) * plotHeight, 0, angstroms, 0, 0, 0 );
			}
			{
				char font[256] = "./FreeMonoBoldOblique.ttf";
				char numb[32];
				int bpos = lc - numletters - 1;
				if ( lc == end - 1 ) bpos = ( ( end - beg ) / numletters ) * numletters + beg - 1;
				sprintf( numb, "%d", offset + bpos );
				image.plot_text( font, 12, 2, posY - plotHeight + 6 /*- ( letter.size() + 1 ) * fontsize - 2*/, 0, numb, 0, 0, 0 );
				bpos = lc;
				sprintf( numb, "%3d", offset + bpos - 1 );
				int poslnum = posX - 32;
				if ( lc == end - 1 ) poslnum += 32;
				poslnum -= strlen( numb ) * 12 - 21;
				image.plot_text( font, 12, poslnum, posY - plotHeight + 6 /*- ( letter.size() + 1 ) * fontsize - 2*/, 0, numb, 0, 0, 0 );
			}
    		if ( lc != end - 1 )
    		{
    			posX  = 0;
    			posY -= interval;
    		}
    	}    	
    	if ( ( lc - beg ) % numletters != numletters - 1 && lc < end - 1 )
    	{
			for ( int vc = 0; vc < value.size(); vc++ )
			{
				double v1 = fmin( value[vc][lc], maxPlot );
				double v2 = fmin( value[vc][ lc + 1 ], maxPlot );
				//if ( value[vc][lc] < maxPlot && value[vc][ lc + 1 ] < maxPlot )
				{
					image.line( posX + width / 2, posY - plotHeight + ( v1 - minPlot ) / ( maxPlot - minPlot ) * plotHeight, 
						posX + width + width / 2, posY - plotHeight + ( v2 - minPlot ) / ( maxPlot - minPlot ) * plotHeight, color[vc][0], color[vc][1], color[vc][2] );
				}
			}
    	}
    	for ( int sc = 0; sc < letter.size(); sc++ )
		{
			char ch[2] = "";
			ch[0] = letter[sc][lc];
			double col[3] = { 0, 0, 0 };
			char font[256] = "./FreeMono.ttf";
			if ( value.size() >= 3 && value[2][lc] - ( value[0][lc] + value[1][lc] ) / 2.2 > threshold )
			{
				col[0] = 0.6;
				strcpy( font, "./FreeMonoBold.ttf" );
			}
			if ( letter.size() >= 2 && letter[0][lc] != letter[1][lc] )
			{
				col[0] = 0;
				strcpy( font, "./FreeMonoBold.ttf" );
			}
			image.plot_text( font, fontsize, posX, posY - plotHeight - fontsize + ( sc + 1 - letter.size() ) * fontsize, 0, ch, col[0], col[1], col[2] );
		}
    	posX += width;
    }
}

int LoadSequence( const char *fname, string& seq )
{ 
	FILE *ifile = fopen( fname, "rt" );
	if ( !ifile ) return 0;
	do
	{
		int c = getc( ifile );
		if ( c == EOF ) break;
		if ( c >= 'A' && c <= 'Z' ) seq.push_back( c );
	}
	while ( 1 );
	fclose( ifile );
	return 1;
}


int LoadData( const char *fname, vector<double>& dat, int format )
{ 
	FILE *ifile = fopen( fname, "rt" );
	if ( !ifile ) return 0;
	char buf[256];
	while ( fgets( buf, 256, ifile ) )
	{
		char *p;
		if ( format == 2 )
		{
			if ( !strtok( buf, " \t\n" ) ) return 0;
			p = strtok( 0, " \t\n" );
		}
		else p = strtok( buf, "\t\n" );
		if ( !p ) return 0;
		dat.push_back( atof( p ) );
	}
	fclose( ifile );
	return 1;
}

int main( int argc, char **argv )
{
	if ( argc < 3 )
	{
		printf( "arguments: {-seq filename} {-dat filename r-g-b} [-beg num] [-end num] [-threshold val] [-offset num] [-maxplot val] [-fmt2] out_png\n" );
		return 1;
	}
	gbPlotPrintout printout;
	string oname;
	printout.threshold = 1.4;
	printout.maxPlot = 10.5;
	int format = 1;
	for ( int ac = 0; ac < argc; ac++ )
	{
		if ( strcmp( argv[ac], "-seq" ) == 0 )
		{
			ac++;
			printout.letter.resize( printout.letter.size() + 1 );
			if ( !LoadSequence( argv[ac], printout.letter.back() ) )
			{
				printf( "can't read sequence\n" );
				return 1;
			}
		}
		else if ( strcmp( argv[ac], "-fmt2" ) == 0 )
		{
			format = 2;
		}
		else if ( strcmp( argv[ac], "-dat" ) == 0 )
		{
			ac++;
			printout.value.resize( printout.value.size() + 1 );
			if ( !LoadData( argv[ac], printout.value.back(), format ) )
			{
				printf( "can't read data %s\n", argv[ac] );
				return 1;
			}
			ac++;
			printout.color.resize( printout.color.size() + 1 );
			char nbuf[256];
			strcpy( nbuf, argv[ac] );
			for ( char *p = strtok( nbuf, "-" ); p; p = strtok( 0, "-" ) ) printout.color.back().push_back( atof( p ) / 256 );
			if ( printout.color.back().size() != 3 )
			{
				printf( "incorrect color\n" );
				return 1;
			}
		}
		else if ( strcmp( argv[ac], "-beg" ) == 0 )
		{
			ac++;
			printout.beg = atoi( argv[ac] );
		}
		else if ( strcmp( argv[ac], "-end" ) == 0 )
		{
			ac++;
			printout.end = atoi( argv[ac] ) + 1;
		}
		else if ( strcmp( argv[ac], "-offset" ) == 0 )
		{
			ac++;
			printout.offset = atoi( argv[ac] );
		}
		else if ( strcmp( argv[ac], "-threshold" ) == 0 )
		{
			ac++;
			printout.threshold = atof( argv[ac] );
		}
		else if ( strcmp( argv[ac], "-maxplot" ) == 0 )
		{
			ac++;
			printout.maxPlot = atof( argv[ac] );
		}
		else oname = argv[ac];
	}
	if ( printout.end > printout.letter[0].size() ) printout.end = printout.letter[0].size();
	pngwriter image(1600, 280 * ( ( printout.end - printout.beg ) / 100 + 1 ) + 1, 1.0, oname.data() );
	printout.DrawPage( image );
	image.close();
	return 0;
}