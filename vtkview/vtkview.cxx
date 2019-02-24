
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <map>
#include <string>
#include <set>

using namespace std;

#include <vtkSmartPointer.h>
#include <vtkVersion.h>

#include <vtkParametricFunctionSource.h>
#include <vtkTupleInterpolator.h>
#include <vtkTubeFilter.h>
#include <vtkParametricSpline.h>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkCamera.h>


static double resolutionScale = 3.;

// Define interaction style
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static KeyPressInteractorStyle* New();
    vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);
 
    virtual void OnKeyPress() 
    {
      // Get the keypress
      vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
 
      // Output the key that was pressed
      std::cout << "Pressed " << key << std::endl;
 
      // Handle an arrow key
      if(key == "Up")
        {
        std::cout << "The up arrow was pressed." << std::endl;
        }
 
      // Handle a "normal" key
      if(key == "a")
        {
		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = 
			vtkSmartPointer<vtkWindowToImageFilter>::New();
		windowToImageFilter->SetInput( rwi->GetRenderWindow() );
		windowToImageFilter->SetMagnification( resolutionScale ); //set the resolution of the output image (3 times the current resolution of vtk render window)
		windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
		//windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
		windowToImageFilter->Update();
		printf( "enter file name to save image [.png]\n" );
	    char fname[256];
	    fgets( fname, 256, stdin );
		if ( strlen( fname ) > 1 )
		{
			strtok( fname, "\n" );
			if ( strlen( fname ) < 4 || strcmp( fname + strlen( fname ) - 4, ".png" ) != 0 ) strcat( fname, ".png" );
			
			vtkSmartPointer<vtkPNGWriter> writer = 
				vtkSmartPointer<vtkPNGWriter>::New();
			writer->SetFileName( fname );
			writer->SetInputConnection(windowToImageFilter->GetOutputPort());
			writer->Write();
			printf( "image saved to %s\n", fname );
		}
		
		rwi->Render();  
		rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Clear();
		rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->ResetCamera();
		rwi->Render(); 
		//std::cout << "The a key was pressed." << std::endl;
        }

     if(key == "c")
        {
		printf( "enter file name to save camera [.cam]\n" );
	    char fname[256];
	    fgets( fname, 256, stdin );
		if ( strlen( fname ) > 1 )
		{
			strtok( fname, "\n" );
			if ( strlen( fname ) < 4 || strcmp( fname + strlen( fname ) - 4, ".cam" ) != 0 ) strcat( fname, ".cam" );
			vtkCamera *cam = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
			FILE *ofile = fopen( fname, "wt" );
			double position[3];
			double focalpoint[3];
			double viewup[3];
			cam->GetPosition( position );
			cam->GetFocalPoint( focalpoint );
			cam->GetViewUp( viewup );
			fprintf( ofile, "%g %g %g\n", position[0], position[1], position[2] );
			fprintf( ofile, "%g %g %g\n", focalpoint[0], focalpoint[1], focalpoint[2] );
			fprintf( ofile, "%g %g %g\n", viewup[0], viewup[1], viewup[2] );
			fclose( ofile );
			printf( "camera attributes saved to %s\n", fname );
		}
        }
     if ( key.size() == 1 && key != "a" && key != "c" && key[0] > 'a' && key[0] <= 'z' )
	 {
		 printf( "a - save image; c - save camera; you pressed %s\n", key.data() );
	 }

      // Forward events
      vtkInteractorStyleTrackballCamera::OnKeyPress();
    }
 
};
vtkStandardNewMacro(KeyPressInteractorStyle);


int main(int argc, char *argv[])
{
	if ( argc == 1 )
	{
		fprintf( stderr, "arguments: pdb corr [cam] [-nomark] [-noscale] [-fragm beg end] [-offset num] [-ignore num] [-title text] [-magnification num=3]\n" );
		return 1;
	}
	bool domark = true;
	bool doscale = true;
	bool camf = false;
	double camattr[3][3];
	char buf[1024];
	int markbeg = 297;
	int markend = 304;
	int moffset = 0;
	int ignorenum = 0;
	string wtitle = argv[1];
	for ( int ac = 1; ac < argc; ac++ )
	{
		if ( strcmp( argv[ac], "-nomark" ) == 0 ) domark = false;
		if ( strcmp( argv[ac], "-noscale" ) == 0 ) doscale = false;
		if ( strcmp( argv[ac], "-fragm" ) == 0 )
		{
			ac++;
			markbeg = atoi( argv[ac] );
			ac++;
			markend = atoi( argv[ac] );
		}
		if ( strcmp( argv[ac], "-offset" ) == 0 )
		{
			ac++;
			moffset = atoi( argv[ac] );
		}
		if ( strcmp( argv[ac], "-ignore" ) == 0 )
		{
			ac++;
			ignorenum = atoi( argv[ac] );
		}
		if ( strcmp( argv[ac], "-magnification" ) == 0 )
		{
			ac++;
			resolutionScale = atof( argv[ac] );
		}
		if ( strcmp( argv[ac], "-title" ) == 0 )
		{
			ac++;
			wtitle = argv[ac];
		}
	}
	if ( argc > 3 && strcmp( argv[ 3 ], "-nomark" ) != 0 && strcmp( argv[3], "-noscale" ) != 0 ) 
	{
		FILE *camfile = fopen( argv[3], "rt" );
		if ( !camfile )
		{
			fprintf( stderr, "can't open camera file %s\n", argv[3] );
			return 1;
		}
		int ac = 0;
		while ( fgets( buf, 1024, camfile ) )
		{
			int cc = 0;
			for ( char *p = strtok( buf, " \n" ); p && cc < 3; p = strtok( 0, " \n" ), cc++ ) camattr[ ac ][ cc ] = atof( p );
			ac++;
		}
		fclose( camfile );
		camf = true;
	}
	FILE *cfile = fopen( argv[2], "rt" );
	
	const int nmut = 7;
	int mutpos[ nmut ] = { 518, 633, 676, 691, 723, 802, 831 };
	int mutflags[ nmut ] = { 2, 1, 1, 1, 1, 2, 2 };
	const int nzn = 2;
	int zncoord[ nzn ][ 3 ] = { { 0 } };
	map<int,int> smutpos;
	for ( int i = 0; i < nmut; i++ ) smutpos[ mutpos[i] ] = i;
	double mutcoord[ nmut ][ 3 ] = { { 0 } };
	
	const int nseg = 35;
	double bestsdist[ nseg ];
	int bestsnum[ nseg ];
	const char *romand[] = { "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x", "xi", "xii", "xiii", "xiv", "xv", "xvi", " xvii", "xviii", "xix", "xx", "xxi", "xxii", "xxiii", "xxiv", "xxv" };
	multimap<double,int> segments;
	if ( !fgets( buf, 1024, cfile ) ) return 1;
	while ( fgets( buf, 1024, cfile ) )
	{
		string s1( strtok( buf, "\t" ) );
		double v = atof( strtok( 0, "\t\n" ) );
		int n = atoi( s1.substr( 2, strcspn( s1.data() + 2, "_" ) ).data() ) / 10;
		segments.insert( pair<double,int>( v, n ) );
	}
	fclose( cfile );
	map<int,double> seln;
	map<int,int> p2rindex0;
	double mincv = segments.begin()->first;
	double maxcv = segments.rbegin()->first;
	if ( doscale )
	{
		mincv = 4e-13;
		maxcv = 9e-12;
	}
	int oc = 0;
	for ( multimap<double,int>::iterator it = segments.begin(); it != segments.end() && oc < nseg; it++, oc++ ) 
	{
		if ( it->second * 10 <= ignorenum ) continue;
		bestsdist[ seln.size() ] = 1000;
		p2rindex0[ it->second ] = seln.size();
		seln[ it->second ] = it->first;
		maxcv = it->first;
		printf( "%d %g\n", it->second, it->first );
	}
		
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> points2 =
		vtkSmartPointer<vtkPoints>::New();
	FILE *ifile = fopen( argv[1], "rt" );
	int pc = moffset;
	int pc2 = 0;
	int pcc = 0;
	map<int,int> p2index;
	map<int,int> p2rindex;
	map<int,double> p2dist;
	double p2dmax = 0;
	double p2dmin = 1000;
	double prev2coord[3] = { 0 };
	int znpos = 0;
	while ( fgets( buf, 1024, ifile ) )
	{
		if ( strncmp( buf, "ATOM", 4 ) == 0 && strncmp( buf + 13, "CA", 2 ) == 0 && !( strncmp( buf + 17, "GTP", 3 ) == 0 || strncmp( buf + 17, "SAH", 3 ) == 0 ) )
		{
			double coord[3];
			coord[0] = atof( buf + 31 );
			coord[1] = atof( buf + 39 );
			coord[2] = atof( buf + 47 );
			int resnum = atoi( buf + 23 );
			if ( resnum < ignorenum ) 
			{
				markbeg--;
				markend--;
				pc++;
				continue;
			}
			points->InsertPoint( pcc, coord[0], coord[1], coord[2] );
			if ( smutpos.find( pc ) != smutpos.end() )
			{
				mutcoord[ smutpos[pc] ][ 0 ] = coord[0];
				mutcoord[ smutpos[pc] ][ 1 ] = coord[1];
				mutcoord[ smutpos[pc] ][ 2 ] = coord[2];
			}
				
			pc++;
			pcc++;
			if ( resnum % 10 == 4 )
			{
				if ( seln.find( resnum / 10 ) != seln.end() )
				{
					points2->InsertPoint( pc2, coord[0], coord[1], coord[2] );
					if ( pc2 > 0 )
					{
						double dist = sqrt( pow( coord[0] - prev2coord[0], 2 ) + pow( coord[1] - prev2coord[1], 2 ) + pow( coord[2] - prev2coord[2], 2 ) );
						p2dmax = fmax( dist, p2dmax );
						p2dmin = fmin( dist, p2dmin );
						p2dist[ pc2 ] = dist;
						printf( "%d %g\n", pc2, dist );
					}
					p2index[ pc2 ] = resnum / 10;
					p2rindex[ p2rindex0[ resnum / 10 ] ] = pc2;
					//printf( "%d %g\n", pc2, seln[ resnum/10] );
					memcpy( prev2coord, coord, 3 * sizeof( double ) );
					pc2++;
				}
			}
		}
		if ( strncmp( buf, "ATOM", 4 ) == 0 && ( strncmp( buf + 12, "ZN", 2 ) == 0 || strncmp( buf + 13, "ZN", 2 ) == 0 ) )
		{
			double coord[3];
			coord[0] = atof( buf + 31 );
			coord[1] = atof( buf + 39 );
			coord[2] = atof( buf + 47 );
			zncoord[ znpos ][ 0 ] = coord[0];
			zncoord[ znpos ][ 1 ] = coord[1];
			zncoord[ znpos ][ 2 ] = coord[2];
			znpos++;
		}
	}
	fclose( ifile );

  // Fit a spline to the points
  vtkSmartPointer<vtkParametricSpline> spline =
    vtkSmartPointer<vtkParametricSpline>::New();
  spline->SetPoints(points);
  vtkSmartPointer<vtkParametricFunctionSource> functionSource =
    vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource->SetParametricFunction(spline);
  functionSource->SetUResolution(10 * points->GetNumberOfPoints());
  functionSource->Update();

  vtkSmartPointer<vtkParametricSpline> spline2 =
    vtkSmartPointer<vtkParametricSpline>::New();
  spline2->SetPoints(points2);
  vtkSmartPointer<vtkParametricFunctionSource> functionSource2 =
    vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource2->SetParametricFunction(spline2);
  functionSource2->SetUResolution(10 * points2->GetNumberOfPoints());
  functionSource2->SetScalarModeToU();
  functionSource2->Update();
/*
  vtkSmartPointer<vtkParametricFunctionSource> functionSourceI2 =
    vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSourceI2->SetParametricFunction(spline2);
  functionSourceI2->SetUResolution(10 * points2->GetNumberOfPoints());
  functionSourceI2->Update();*/
  // Interpolate the scalars
  double rad;
  vtkSmartPointer<vtkTupleInterpolator> interpolatedRadius =
    vtkSmartPointer<vtkTupleInterpolator> ::New();
  interpolatedRadius->SetInterpolationTypeToLinear();
  interpolatedRadius->SetNumberOfComponents(1);
  for ( int i = 0; i < points->GetNumberOfPoints(); i++ )
  {
	rad = .1; 
	if ( i >= markbeg && i <= markend ) rad = .4;
	interpolatedRadius->AddTuple(i,&rad);
  }
  vtkSmartPointer<vtkTupleInterpolator> interpolatedRadius2 =
    vtkSmartPointer<vtkTupleInterpolator> ::New();
  interpolatedRadius2->SetInterpolationTypeToLinear();
  interpolatedRadius2->SetNumberOfComponents(1);
  for ( int i = 0; i < points2->GetNumberOfPoints(); i++ )
  {
	rad = 1.2; interpolatedRadius2->AddTuple(i,&rad);
  }
  /*
  vtkSmartPointer<vtkTupleInterpolator> interpolatedIndex2 =
    vtkSmartPointer<vtkTupleInterpolator> ::New();
  interpolatedIndex2->SetInterpolationTypeToLinear();
  interpolatedIndex2->SetNumberOfComponents(1);
  for ( int i = 0; i < points2->GetNumberOfPoints(); i++ )
  {
	rad = i; interpolatedIndex2->AddTuple(i,&rad);
  }
  */
	/*
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint(0,1,0,0);
  points->InsertPoint(1,2,0,0);
  points->InsertPoint(2,3,1,0);
  points->InsertPoint(3,4,1,0);
  points->InsertPoint(4,5,0,0);
  points->InsertPoint(5,6,0,0);

  // Fit a spline to the points
  vtkSmartPointer<vtkParametricSpline> spline =
    vtkSmartPointer<vtkParametricSpline>::New();
  spline->SetPoints(points);
  vtkSmartPointer<vtkParametricFunctionSource> functionSource =
    vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource->SetParametricFunction(spline);
  functionSource->SetUResolution(10 * points->GetNumberOfPoints());
  functionSource->Update();

  // Interpolate the scalars
  double rad;
  vtkSmartPointer<vtkTupleInterpolator> interpolatedRadius =
    vtkSmartPointer<vtkTupleInterpolator> ::New();
  interpolatedRadius->SetInterpolationTypeToLinear();
  interpolatedRadius->SetNumberOfComponents(1);
  rad = .2; interpolatedRadius->AddTuple(0,&rad);
  rad = .2; interpolatedRadius->AddTuple(1,&rad);
  rad = .2; interpolatedRadius->AddTuple(2,&rad);
  rad = .1; interpolatedRadius->AddTuple(3,&rad);
  rad = .1; interpolatedRadius->AddTuple(4,&rad);
  rad = .1; interpolatedRadius->AddTuple(5,&rad);
	*/
  // Generate the radius scalars
  vtkSmartPointer<vtkDoubleArray> tubeRadius =
    vtkSmartPointer<vtkDoubleArray>::New();
  unsigned int n = functionSource->GetOutput()->GetNumberOfPoints();
  unsigned int n1p = points->GetNumberOfPoints();
  unsigned int n2p = points2->GetNumberOfPoints();
  tubeRadius->SetNumberOfTuples(n);
  tubeRadius->SetName("TubeRadius");
  double tMin = interpolatedRadius->GetMinimumT();
  double tMax = interpolatedRadius->GetMaximumT();
  double r;
  for (unsigned int i = 0; i < n; ++i)
    {
 	int res = ( double( i ) * ( n1p - 1 ) ) / n; 
	int rf = 255, gf = 255, bf = 255;
	int flag = ( res >= markbeg && res <= markend );
	if ( flag ) r = 0.3;
	else r = 0.2;
	tubeRadius->SetTuple1(i, r);
    }

  vtkSmartPointer<vtkDoubleArray> tubeRadius2 =
    vtkSmartPointer<vtkDoubleArray>::New();
  unsigned int n2 = functionSource2->GetOutput()->GetNumberOfPoints();
  printf( "n2 = %u, n = %u, n1p = %u, n2p = %u markbeg = %d markend = %d\n", n2, n, n1p, n2p, markbeg, markend );
  tubeRadius2->SetNumberOfTuples(n2);
  tubeRadius2->SetName("TubeRadius2");
  double tMin2 = interpolatedRadius2->GetMinimumT();
  double tMax2 = interpolatedRadius2->GetMaximumT();
  for (unsigned int i = 0; i < n2; ++i)
    {
    //double t = (tMax2 - tMin2) / (n2 - 1) * i + tMin2;
    //interpolatedRadius2->InterpolateTuple(t, &r);
		double u2 = double( i ) / n2;
		double su[3];
		double sp[3];
		double sd[9];
		su[0] = u2;
		spline2->Evaluate( su, sp, sd );
		for ( int k = 0; k < nseg; k++ )
		{
			const double *coord = points2->GetPoint( k );
			double d = 0;
			for ( int j = 0; j < 3; j++ ) d += ( coord[j] - sp[j] ) * ( coord[j] - sp[j] );
			d = sqrt( d );
			if ( bestsdist[k] > d )
			{
				bestsdist[k] = d;
				bestsnum[k] = i;
			}
		}
	}
			
  for (unsigned int i = 0; i < n2; ++i)
    {
		int cj = 0;
		for ( int j = 0; j < n2p - 1; j++ )
		{
			if ( bestsnum[j] <= i && bestsnum[ j + 1 ] >= i )
			{
				cj = j;
				break;
			}
		}
	double r0 = 1.2;
	double dist = p2dist[ cj + 1 ];
	double ir1 = 5 * double( i - bestsnum[cj] ) / ( bestsnum[ cj + 1 ] - bestsnum[cj] );
	double ir2 = 5 * double( bestsnum[ cj + 1 ] - i ) / ( bestsnum[ cj + 1 ] - bestsnum[cj] );
	double rc = fmax( exp( - ir1 * ir1 ), exp( -ir2 * ir2 ) );
	double r = r0 * ( 1. - ( 1. - rc ) * dist / p2dmax );
	printf( "%d %d %g %g %g %g\n", i, cj, dist, ir1, ir2, r );
	if ( r < 0.3 ) r = 0.3;
	//r = 1.2;
    tubeRadius2->SetTuple1(i, r);
    }
/*
  vtkSmartPointer<vtkDoubleArray> tubeIndex2 =
    vtkSmartPointer<vtkDoubleArray>::New();
  unsigned int nI2 = functionSourceI2->GetOutput()->GetNumberOfPoints();
  tubeIndex2->SetNumberOfTuples(nI2);
  tubeIndex2->SetName("TubeIndex2");
  double tIMin2 = interpolatedIndex2->GetMinimumT();
  double tIMax2 = interpolatedIndex2->GetMaximumT();
  
  for (unsigned int i = 0; i < nI2; ++i)
    {
    double t = (tIMax2 - tIMin2) / (nI2 - 1) * i + tIMin2;
    interpolatedIndex2->InterpolateTuple(t, &r);
    tubeIndex2->SetTuple1(i, r);
    }
    */
  // Add the scalars to the polydata
  vtkSmartPointer<vtkPolyData> tubePolyData =
    vtkSmartPointer<vtkPolyData>::New();
  tubePolyData = functionSource->GetOutput();
  tubePolyData->GetPointData()->AddArray(tubeRadius);
  tubePolyData->GetPointData()->SetActiveScalars("TubeRadius");

  vtkSmartPointer<vtkUnsignedCharArray> rcolors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  rcolors->SetName("RColors");
  rcolors->SetNumberOfComponents(3);
  rcolors->SetNumberOfTuples(n);
  for (int i = 0; i < n ;i++)
  {
	  //double pr = ( seln[ p2index[ i / 10 ] ] - mincv ) / maxcv - mincv;
	  //printf( "%g\n", seln[ p2index[ i / 10 ] ] );
	  //double pr = ( log( seln[ p2index[ i / 10 ] ] ) - log( mincv ) ) / ( log( maxcv ) - log( mincv )  );
	  //double pr = ( log( seln[ p2index[ i / 10 ] ] ) ) / ( log( maxcv ) );
	int res = ( double( i ) * ( n1p - 1 ) ) / n; 
	int rf = 255, gf = 255, bf = 255;
	int flag = ( res >= markbeg && res <= markend );
	if ( flag )
	{
		rf = 255;
		gf = 10;
		bf = 10;
	}
	if ( seln.find( res / 10 ) != seln.end() && flag == 255 && domark && false ) 
	{
		rf = 200;
		gf = 10;
		bf = 10;
	}
    rcolors->InsertTuple3(i,
                       rf,
                       gf,
                       bf );
  }
  tubePolyData->GetPointData()->AddArray(rcolors);

  // Create the tubes
  vtkSmartPointer<vtkPolyData> tubePolyData2 =
    vtkSmartPointer<vtkPolyData>::New();
  tubePolyData2 = functionSource2->GetOutput();
  tubePolyData2->GetPointData()->AddArray(tubeRadius2);
  tubePolyData2->GetPointData()->SetActiveScalars("TubeRadius2");

  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetName("Colors");
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples(n2);
  for (int i = 0; i < n2 ;i++)
  {
	  //double pr = ( seln[ p2index[ i / 10 ] ] - mincv ) / maxcv - mincv;
	  //printf( "%g\n", seln[ p2index[ i / 10 ] ] );
		int cj = 0;
		for ( int j = 0; j < n2p - 1; j++ )
		{
			if ( bestsnum[j] <= i && bestsnum[ j + 1 ] >= i )
			{
				cj = j;
				break;
			}
		}
		double ir1 = double( i - bestsnum[cj] ) / ( bestsnum[ cj + 1 ] - bestsnum[cj] );
		double ir2 = double( bestsnum[ cj + 1 ] - i ) / ( bestsnum[ cj + 1 ] - bestsnum[cj] );
		if ( ir1 < ir2 ) { ir1 = 0; ir2 = 1; }
		else { ir1 = 1; ir2 = 0; }
		double cval = ( 1 - ir1 ) * seln[ p2index[ cj ] ] + ( 1 - ir2 ) * seln[ p2index[ cj + 1 ] ];
		
	  double pr = ( log( cval ) - log( mincv ) ) / ( log( maxcv ) - log( mincv )  );
	  //double pr = ( log( seln[ p2index[ i / 10 ] ] ) ) / ( log( maxcv ) );
    colors->InsertTuple3(i,
                       0,
                       int(255. * pr ),
                       int(255. * (1 - pr) ) );
  }
  tubePolyData2->GetPointData()->AddArray(colors);
  
  /*
  vtkSmartPointer<vtkPolyData> tubePolyDataI2 =
    vtkSmartPointer<vtkPolyData>::New();
  tubePolyDataI2 = functionSourceI2->GetOutput();
  tubePolyDataI2->GetPointData()->AddArray(tubeIndex2);
  tubePolyDataI2->GetPointData()->SetActiveScalars("TubeIndex2");
  */
  vtkSmartPointer<vtkTubeFilter> tuber =
    vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  tuber->SetInput(tubePolyData);
#else
  tuber->SetInputData(tubePolyData);
#endif
  tuber->SetNumberOfSides(20);
  tuber->SetVaryRadiusToVaryRadiusByAbsoluteScalar();

  vtkSmartPointer<vtkTubeFilter> tuber2 =
    vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  tuber2->SetInput(tubePolyData2);
#else
  tuber2->SetInputData(tubePolyData2);
#endif
  tuber2->SetNumberOfSides(20);
  tuber2->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
  /*
  //tubePolyData2->GetPointData()->SetActiveScalars("TubeIndex2");
  //--------------
  // Setup actors and mappers
  vtkSmartPointer<vtkPolyDataMapper> lineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  lineMapper->SetInput(tubePolyData);
#else
  lineMapper->SetInputData(tubePolyData);
#endif
  lineMapper->SetScalarRange(tubePolyData->GetScalarRange());

  vtkSmartPointer<vtkPolyDataMapper> lineMapper2 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  lineMapper2->SetInput(tubePolyDataI2);
#else
  lineMapper2->SetInputData(tubePolyDataI2);
#endif
	lineMapper2->SetScalarRange(tubePolyDataI2->GetScalarRange());
	*/
  vtkSmartPointer<vtkPolyDataMapper> tubeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  tubeMapper->SetInputConnection(tuber->GetOutputPort());
  tubeMapper->SetScalarRange(tubePolyData->GetScalarRange());
  tubeMapper->ScalarVisibilityOn();
  tubeMapper->SetScalarModeToUsePointFieldData();
  tubeMapper->SelectColorArray("RColors");

  vtkSmartPointer<vtkPolyDataMapper> tubeMapper2 =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  tubeMapper2->SetInputConnection(tuber2->GetOutputPort());
  tubeMapper2->SetScalarRange(tubePolyData2->GetScalarRange());
  tubeMapper2->ScalarVisibilityOn();
  tubeMapper2->SetScalarModeToUsePointFieldData();
  tubeMapper2->SelectColorArray("Colors");

  //vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
  //lineActor->SetMapper(lineMapper);
  vtkSmartPointer<vtkActor> tubeActor = vtkSmartPointer<vtkActor>::New();
  tubeActor->SetMapper(tubeMapper);

  //vtkSmartPointer<vtkActor> lineActor2 = vtkSmartPointer<vtkActor>::New();
  //lineActor2->SetMapper(lineMapper2);
  vtkSmartPointer<vtkActor> tubeActor2 = vtkSmartPointer<vtkActor>::New();
  tubeActor2->SetMapper(tubeMapper2);

  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetWindowName( wtitle.data() );
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  
  vtkSmartPointer<KeyPressInteractorStyle> style = 
    vtkSmartPointer<KeyPressInteractorStyle>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  style->SetCurrentRenderer(renderer);
  
  //renderer->AddActor(lineActor);
  renderer->AddActor(tubeActor);
  //renderer->AddActor(lineActor2);
  renderer->AddActor(tubeActor2);
  
  for (int mc = 0; mc < nmut; mc++ )
  {
		/*
		vtkSmartPointer<vtkUnsignedCharArray> scolors =
			vtkSmartPointer<vtkUnsignedCharArray>::New();
		scolors->SetName("SColors");
		scolors->SetNumberOfComponents(3);
		scolors->SetNumberOfTuples(1);
		scolors.InsertNextTuple3(255,0,0);
		*/
		vtkSmartPointer<vtkSphereSource> sphereSource = 
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter( mutcoord[mc][0], mutcoord[mc][1], mutcoord[mc][2] );
		sphereSource->SetRadius(1.7);
		/*
		sp = sphereSource->GetOutput();
		sphereSource->Update();
		sp->GetCellData()->SetScalars( scolors );
		appendData.AddInputData(sp )
		*/
		
		vtkSmartPointer<vtkPolyDataMapper> smapper = 
			vtkSmartPointer<vtkPolyDataMapper>::New();
		smapper->SetInputConnection(sphereSource->GetOutputPort());
		//smapper->SetColorModeToDirectScalars();
		//smapper->SelectColorArray("SColors");

		vtkSmartPointer<vtkActor> sactor = 
		vtkSmartPointer<vtkActor>::New();
		sactor->SetMapper(smapper);
		if ( mutflags[ mc ] == 1 )
		{
			sactor->GetProperty()->SetColor(0.1, 0.1, 0.1);
		}
		else
		{
			sactor->GetProperty()->SetColor(0.1, 0.1, 0.1);
		}
			
		renderer->AddActor(sactor);
  }
  for (int zc = 0; zc < nzn; zc++ )
  {
		/*
		vtkSmartPointer<vtkUnsignedCharArray> scolors =
			vtkSmartPointer<vtkUnsignedCharArray>::New();
		scolors->SetName("SColors");
		scolors->SetNumberOfComponents(3);
		scolors->SetNumberOfTuples(1);
		scolors.InsertNextTuple3(255,0,0);
		*/
		vtkSmartPointer<vtkSphereSource> sphereSource = 
			vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter( zncoord[zc][0], zncoord[zc][1], zncoord[zc][2] );
		sphereSource->SetRadius(2.7);
		/*
		sp = sphereSource->GetOutput();
		sphereSource->Update();
		sp->GetCellData()->SetScalars( scolors );
		appendData.AddInputData(sp )
		*/
		
		vtkSmartPointer<vtkPolyDataMapper> smapper = 
			vtkSmartPointer<vtkPolyDataMapper>::New();
		smapper->SetInputConnection(sphereSource->GetOutputPort());
		//smapper->SetColorModeToDirectScalars();
		//smapper->SelectColorArray("SColors");

		vtkSmartPointer<vtkActor> sactor = 
		vtkSmartPointer<vtkActor>::New();
		sactor->SetMapper(smapper);
		sactor->GetProperty()->SetColor(0.57, 0, 0.82);
		renderer->AddActor(sactor);
  }
  if ( domark ) for (int sc = 0; sc < points2->GetNumberOfPoints(); sc++ )
  {
	  const double *coord = points2->GetPoint( p2rindex[ sc ] );
	  double ccoord[3];
	  ccoord[0] = coord[0];
	  ccoord[1] = coord[1];
	  ccoord[2] = coord[2];
	  const char *crdigit = romand[ sc ];
	  ccoord[0] += 1.6;
	  for ( int sc1 = 0; sc1 < strlen( crdigit ); sc1++ )
	  {
			char cd = crdigit[ sc1 ];
			vtkSmartPointer<vtkSphereSource> sphereSource = 
				vtkSmartPointer<vtkSphereSource>::New();
			sphereSource->SetCenter( ccoord[0], ccoord[1], ccoord[2] );
			sphereSource->SetRadius( cd == 'x' ? 1.2 : ( ( cd == 'v' ) ? 0.9 : 0.6 ) );
			vtkSmartPointer<vtkPolyDataMapper> smapper = 
				vtkSmartPointer<vtkPolyDataMapper>::New();
			smapper->SetInputConnection(sphereSource->GetOutputPort());
			vtkSmartPointer<vtkActor> sactor = 
			vtkSmartPointer<vtkActor>::New();
			sactor->SetMapper(smapper);
			sactor->GetProperty()->SetColor( 0.1, 0.1, (cd == 'v' ) ? 0.8: 0.1 );
			renderer->AddActor(sactor);
			ccoord[0] += 1.6;
	  }
  }
  
  //renderer->SetBackground(0.78, 0.78, 0.78);
  renderer->SetBackground(1., 1., 1.);
  if ( camf )
  {
	  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
	  for ( int i = 0; i < 3; i++ ) for( int k = 0; k < 3; k++ ) printf( "%g ", camattr[i][k] );
	  camera->SetPosition( camattr[0] );
	  camera->SetFocalPoint( camattr[1] );
	  camera->SetViewUp( camattr[2] );
	  printf( "\ncamera configured\n" );
	  renderer->SetActiveCamera(camera);
  }
  renderWindow->Render();
  
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
