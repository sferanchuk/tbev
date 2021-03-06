#!/usr/bin/python

# Turn on debug mode.
import sys
tfp = '/home/sferanchuk/.local/lib/python2.7/site-packages'
if tfp in sys.path:
	sys.path.remove( tfp )
import cgi
import cgitb
cgitb.enable()
import csv
import numpy as np
#from sklearn import manifold
#from sklearn import decomposition
#import scipy.stats as stats
#from scipy.spatial import distance
import json
import os
import collections
import operator
import math
from subprocess import Popen, PIPE
import sqlite3

print "Content-Type: text/html\r\n\r\n" 
print """
<head>
        <title>Bublle Chart</title>
<style>

circle {
  fill: rgb(31, 119, 180);
  fill-opacity: 0;
  stroke: rgb(31, 119, 180);
  stroke-width: 1px;
}

.leaf circle {
  fill: #ff7f0e;
  fill-opacity: 1;
}

text {
  font: 10px sans-serif;
  text-anchor: middle;
}

.node {
  font: 10px sans-serif;
  line-height: 12px;
  overflow: hidden;
  position: absolute;
  text-indent: 2px;
}

</style>
	  <script type="text/javascript" src="js/d3.v4.min.js"></script>
	<script type="text/javascript" src=js/d3-scale-chromatic.v0.3.min.js></script>
	<script type="text/javascript" src="js/FileSaver.js"></script>
 	  <script type="text/javascript" src="js/html2canvas.js"></script>
      <script type="text/javascript" src="js/html2canvas.svg.js"></script>
	  <script type="text/javascript" src="js/canvas-to-blob.min.js"></script>
</head>
<body>
"""


dir1 = "/home/sferanchuk/disk2/ns5/dal/segments"
dir2 = "/home/sferanchuk/disk2/ns5/d_dal/segments"

ns5dal = "GGSEGDTLGDLWKRKLNGCTKEEFFAYRRTGILETERDKARELLRRGETNMGLAVSRGTAKLAWLEERGYATLKGEVVDLGCGRGGWSYYAASRPAVMSVKAYTIGGKGHETPKMVTSLGWNLIKFRAGMDVFSMQPHRADTIMCDIGESNPDAVVEGERTRKVILLMEQWKNRNPTATCVFKVLAPYRPEVIEALHRFQLQWGGGLVRTPFSRNSTHEMYYSTAVTGNIVNSVNIQSRKLLARFGDQRGPTRVPELDLGVGTRCVVLAEDKVKEKDVQERISALREQYGETWHMDREHPYRTWQYWGSYRTAPTGSAASLINGVVKLLSWPWNAREDVVRMAMTDTTAFGQQRVFKEKVDTKAQEPQPGTKVIMRAVNDWILERLARKSKPRMCSREEFIAKVKSNAALGAWSDEQNRWSSAKEAVEDPAFWQLVDEERERHLAGRCAHCVYNMMGKREKKLGEFGVAKGSRAIWYMWLGSRFLEFEALGFLNEDHWASRGSSGSGVEGISLNYLGWHLKGLSTLEGGLFYADDTAGWDTKVTNADLEDEEQLLRYMEGEHKQLAATIMQKAYHAKVVKVARPSRDGGCIMDVITRRDQRGSGQVVTYALNTLTNIKVQLIRMMEGEGVIEASDAHNPRLLRVERWLRDHGEERLGRMLVSGDDCVVRPVDDRFSRALYFLNDMAKTRKDIGEWEHSVGFSNWEEVPFCSHHFHELVMKDGRALIVPCRDQDELVGRARVSPGCGWSVRETACLSKAYGQMWLLSYFHRRDLRTLGLAICSAVPVDWVPTGRTTWSIHASGSWMTTEDMLDVWNRVWILDNPFMHSKEKIAEWRDVPYLPKSHDMLCSSLVGRKERAEWAKNIWGAVEKVRKMIGQEKFKDYLSCMDRHDLHWELKLESSII"
models = [ "hpol", "lpol", "hpol2", "lpol2" ]
smodels = [ "dal", "pri", "d_dal", "d_pri" ]
mfiles = [ "dal-znmg-dyn3-f.pdb", "pri-znmg-430-480-f.pdb", "ns5_d_dalc_dyn3_f.pdb",  "ns5_d_pric_dyn2_f.pdb"  ]
avgmodels = [ "all", "mt", "int", "rdrp" ]
segments = [ "segmlist1", "segmlist1", "segmlist2", "segmlist2", "segmlist2" ]
snames = { "321-331" : "bN", "331-341" : "bN1", "351-361" : "bN2", "361-371" : "aN", "371-381" : "aN1", "381-391" : "aN2", "391-401" : "aN3", "271-281" : "Lin", "261-271" : "Piv","311-321" : "L1",  "341-351" : "L2", "411-421" : "L3", "511-521" : "1s", "631-641" : "2s", "671-681" : "3s", "691-701" : "4s", "721-731" : "5s", "801-811" : "6s", "831-841" : "7s", "101-111" : "SAH", "291-301" : "R", "301-311" : "R1", "51-61" : "1", "141-151": "2", "181-191" : "3", "211-221" : "4",  "241-251": "FGDQ", "401-411" : "G", "41-51": "LRG", "441-451": "Zn1", "451-461": "F1", "461-471": "F2", "471-481" : "F3", "481-491" : "F4", "531-541" : "A", "601-611" : "B", "661-671": "C", "691-701": "D", "711-721": "Zn2", "841-851" : "Zn2", "721-731": "E", "751-761" : "H", "781-791": "Pr1", "791-801": "Pr2", "801-811": "Pr3" }
groups = [ [ [ "hpol" ], [ "lpol" ], [ "hpol2" ], [ "lpol2" ] ], [ [ "hpol", "lpol" ], [ "hpol2", "lpol2" ] ], [ [ "hpol", "hpol2" ], [ "lpol", "lpol2" ] ] ]
groupnames = [ "all models", "state0-state1", "h-l"  ]
gcnames = [ [ "Hpol", "Lpol", "Hpol2", "Lpol2" ], [ "State 1", "State 2" ], [ "H", "L" ] ]
maxitems = 100
gnum = 0
dglabels = "yes"


form = cgi.FieldStorage()
sfragment = form.getvalue( "sfragment" )
distr = form.getvalue( "distr", "zipf" )
nmode = int( form.getvalue( "nmode", "1" ) )
gnum = int( form.getvalue( "gnum", "0" ) )
dglabels = form.getvalue( "dglabels", "yes" )
bnames = form.getvalue( "bnames", "yes" )
truncate = int( form.getvalue( "truncate", "5" ) )
maxitems = int( form.getvalue( "maxitems", "100" ) )
cscheme = form.getvalue( "cscheme", "schemeDark2" )
numb = form.getvalue( "numb", "from 11" )
mchecked = [ 0 ] * 5
for mc in range( len( models ) ):
	if form.getvalue( models[mc], "0" ) != "0":
		mchecked[ mc ] = 1
#dgroup = form.getvalue( "dgroup", "none" )
#dtlabels = form.getvalue( "dtlabels", "no" )
#order2 = form.getvalue( "order2", "none" )
#labels = form.getvalue( "labels", "none" )
#dfilter = form.getvalue( "filter", "none" )
#level = form.getvalue( "level" )
#dnorm = form.getvalue( "dnorm", "percent" )
#dprestype = form.getvalue( "dprestype", "bubble" )
dbname = "nmabubbles.db"
if not os.path.isfile( dbname ):
	conn = sqlite3.connect( dbname )
	cur = conn.cursor()
	cur.execute( "create table segments( nmode text, model text, output text);" )
	conn.commit()
	conn.close()



if True:
	
	nedata = {}
	sgnames = {}
	sbnames = {}
	#print "%d %d<br>" % ( beg, beg + len( sfragment ) )
	for mc in range( len( models ) ):
		if mchecked[ mc ] == 0:
			continue
		moffset = 0
		if mc < 2:
			cdir = dir1
		else:
			cdir = dir2
			moffset = 4
		if numb == "from 11":
			moffset = -moffset
			
		conn = sqlite3.connect( dbname )
		cur = conn.cursor()
		cur.execute( "select output from segments where nmode = ? and model = ?", ( str( nmode ) + numb, models[ mc ] ) )
		dbres = cur.fetchall()
		conn.close()
		if len( dbres ) > 0:
			output = dbres[0][0]
		else:
			output = Popen( [ '/mnt/disk2/www/mdtraj/nma.sh', mfiles[ mc ], str( nmode + 6 ), str( moffset ) ], stdout=PIPE).communicate()[0]
			conn = sqlite3.connect( dbname )
			cur = conn.cursor()
			cur.execute( "INSERT INTO segments VALUES ('%s','%s', '%s')" % ( str( nmode ) + numb, models[mc], output ) )
			conn.commit()
			conn.close()
		cedata = {}
		sres = output.split( "\n" )
		cedata = {}
		for sline in sres:
			sfields = sline.split( "\t" )
			if len( sfields ) < 2:
				continue
			sn = sfields[0].split( "_" )
			if len( sn ) < 2:
				continue
			#if int( sn[1] ) == 1 or int( sn[1] ) > 880:
			#	continue
			rpos = int( sn[1] )
			if not sfields[0] in sgnames:
				if rpos < 237:
					sgnames[ sfields[ 0 ] ] = "Mt"
				elif rpos < 310:
					sgnames[ sfields[ 0 ] ] = "Int"
				elif rpos <= 495 or ( rpos > 544 and rpos < 609 ):
					sgnames[ sfields[ 0 ] ] = "Fing"
				elif rpos < 712:
					sgnames[ sfields[ 0 ] ] = "Palm"
				else:
					sgnames[ sfields[ 0 ] ] = "Thumb"
				sname = sn[1] + "-" + sn[2] 
				if bnames == "yes" and sname in snames:
					sname = snames[ sname ]
				sbnames[ sfields[0] ] = sname
			cedata[ sfields[ 0 ] ] = float( sfields[ 1 ] )
		sedata = collections.OrderedDict( sorted( cedata.items(), key=operator.itemgetter(1), reverse=False ) )
		oc = 1
		for ik, ij in sedata.iteritems():
			if oc <= maxitems:
				if not ik in nedata:
					nedata[ik] = [ 0 ] * len( models )
				v = 0
				if distr == "zipf":
					v = 1. / oc
				elif distr == "pow1.5":
					v = math.pow( oc, -1.5 )
				elif distr == "lognormal0.5":
					v = math.exp( -0.5 * math.log( oc ) * math.log( oc ) )
				elif distr == "natural":
					v = 1000 / ij
				nedata[ ik ][ mc ] = v
				oc += 1
			
	edata = {}
	for key in nedata:
		crow = [ 0 ] * len( groups[ gnum ] )
		for i in range( len( groups[ gnum ] ) ):
			for gkey in groups[ gnum ][ i ]:
				crow[ i ] += nedata[ key ][ models.index( gkey ) ] / len( groups[ gnum ][ i ] )
		edata[ key ] = crow
	#print edata
	#aedata = np.array( edata, dtype=float )
	#aenorm = np.sum( aedata, axis=1 )
	#aedata /= aenorm.reshape( len(edata), 1 )

print "<form method=\"get\" action=\"nmabubbles.py\">"
if True:
	selstr = '<select name="{0}">\n{1}</select>\n'
	optstr = '<option value="{0}" {1} >{0}</option>\n'
	optstr2 = '<option value="{0}" {1} >{2}</option>\n'
	print "<br>Normal mode "
	print selstr.format( "nmode", ''.join( optstr.format(v, "selected" if v == str( nmode ) else "" ) for v in [ "1", "2", "3", "4", "5", "6" ] ) )
	print "<br>Group "
	print selstr.format( "gnum", ''.join( optstr2.format( str( v ), "selected" if v == gnum else "", groupnames[ v ] ) for v in range( len( groups ) ) ) )
	print "<br>Truncate "	
	print selstr.format( "truncate", ''.join( optstr.format( str( v ), "selected" if v == truncate else "" ) for v in [ 0, 3, 5, 10, 20 ] ) )
	print "<br>Distr "	
	print selstr.format( "distr", ''.join( optstr.format(v, "selected" if v == distr else "" ) for v in [ "zipf", "pow1.5", "lognormal0.5", "natural" ] ) )
	print "<br>Group labels "
	print selstr.format( "dglabels", ''.join( optstr.format(v, "selected" if v == dglabels else "" ) for v in [ "no", "yes" ] ) )
	print "<br>Max items "
	print selstr.format( "maxitems", ''.join( optstr.format(v, "selected" if v == str( maxitems ) else "" ) for v in [ "20", "25", "30", "35", "45", "60", "100" ] ) )
	print "<br>bnames "
	print selstr.format( "bnames", ''.join( optstr.format(v, "selected" if v == bnames else "" ) for v in [ "yes", "no" ] ) )
	print "<br>numbering "
	print selstr.format( "numb", ''.join( optstr.format(v, "selected" if v == numb else "" ) for v in [ "from 11", "from exact" ] ) )
	print "<br>cscheme "
	print selstr.format( "cscheme", ''.join( optstr.format(v, "selected" if v == cscheme else "" ) for v in [ "schemeCategory10", "schemeAccent", "schemeDark2", "schemePastel2", "schemeSet2", "schemeSet1", "schemePastel1", "schemePaired" ] ) )
	print "<p>"
	for mc in range( len( models ) ):
		cs=""
		if mchecked[mc] == 1:
			cs = "checked"
		print "%s <input type=checkbox name=\"%s\" value=\"1\" %s><br>" % ( gcnames[ 0 ][ mc ], models[mc], cs )

print """
<br><input type=submit name="ok" value="ok">
</form>
"""
if len( edata ) == 0:
	print "</body></html>"
	sys.exit( 0 )
print """
<button id="generate">Save as SVG</button><br>
<button id="gpng">Save as PNG</button><br>
<div class='content'>
          <!-- /the chart goes here -->
</div>
"""
print "<svg width=960 height=960 id=normal></svg>"
print "<svg width=2400 height=2400 id=highres></svg>"
print "<script type=\"text/javascript\">"
print "var height = 600;"
print "var outfn = \"mdbubble_%s_%d\";" % ( sfragment, gnum ) 
print "var tfont = \"12px sans-serif\";"
#print "var ilevel = " + str( ilevel ) + ";"
print "var dglabels = \"" + dglabels + "\";"
#print "var maxdata = " + str( int( aenorm ) ) + ";"
#print "var nsamples = " + str( len( edata) ) + ";"
#print "var maxidsize = " + str( maxidsize ) + ";"
print "var root =  { name: \"root\", children: ["
gsize = len( groups[ gnum ] )
gscnt = 0
ksize = len( edata )
for gkey in range( len( groups[ gnum ] ) ):
	gn = ""
	if gnum != 0 or mchecked[ models.index( groups[ gnum ][ gkey ][ 0 ] ) ] == 1:
		gn = gcnames[ gnum ][ gkey ]
	print "{ name: \"%s\", children: [ " % gn 
	kcnt = 0
	for ckey in edata:
		csize = int( edata[ ckey ][ gkey ] * 500 )
		#if dnorm == "percent":
		#	csize = int( 10000 * aedata[ gscnt ][ kdict[ ckey ] ] )
		#pdig = [ch.isdigit() for ch in dname].index( True )
		print "{ name: \"%s\", order: \"%s\", size: %d } " % ( sbnames[ ckey ], sgnames[ ckey ], csize )
		if kcnt + 1 < ksize:
			print ","
		kcnt = kcnt + 1
	print "] }"
	if gscnt + 1 < gsize:
		print ","
	gscnt = gscnt + 1
print "] };"
print "var color = d3.scaleOrdinal().domain( [ 'Mt', 'Int', 'Fing', 'Palm', 'Thumb' ] ).range( [\"#63cccd\", \"#fad502\", \"#345a7a\", \"#37c437\", \"#cf9e9e\" ] );" 
print """

// ************** Generate the tree diagram	 *****************
var margin = {top: 100, right: 15, bottom: 100, left: 60}
      , width = 1200 - margin.left - margin.right
     , height = 500 - margin.top - margin.bottom;
//var margin = {top: 20, right: 120, bottom: 20, left: 120},
//	width = 960 - margin.right - margin.left,
//	height = 500 - margin.top - margin.bottom;


root = d3.hierarchy(root)
      .sum(function(d) { return d.size; })
      .sort(function(a, b) { return b.value - a.value; });


function draw( val_id ) {


var svg = d3.select( val_id ),
    //.attr('id', "mainsvg" ),
    diameter = +svg.attr("width");
var fs = diameter / 120;

   var g = svg.append("g").attr("transform", "translate(2,2)"),
    format = d3.format(",d");

  

  var pack = d3.pack()
    .size([diameter - 4, diameter - 4]);




  var node = g.selectAll(".node")
    .data(pack(root).descendants())
    .enter().append("g")
      .attr("class", function(d) { return d.children ? "node" : "leaf node"; })
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

  node.append("title")
	  .style("font", tfont )
	  .style( "text-anchor", "middle" )
      .text(function(d) { return d.data.name; });

  node.append("circle")
      .attr("r", function(d) { return d.r; })
      .style( "stroke", "rgb(31, 119, 180)" )
	  .style( "stroke-width", function(d) { 
			if ( d.data.name == "root" ) 
			{
				return "0px";
			}
			else if ( d.children )
			{
				return "1px";
			}
			else return "1px"; 
			})
      .style( "fill-opacity", function(d) { 
			if ( d.data.name == "root" ) 
			{
				return "0.0";
			}
			else if ( d.children )
			{
				return "1";
			}
			else return "1"; 
			})
      .style("fill", function(d) { 
			if ( d.data.name == "root" ) 
			{
				return d3.hsl( 0.3, 0.3, 0.3 );
			}
			else if ( d.children )
			{
				return d3.hsl( 0, 0, 0.4 );
			}
			else return color(d.data.order); 
			});

  node.filter(function(d) { return !d.children; }).append("text")
      .attr("dy", "0.3em")
      .style("font", fs * 1.3 + "px sans-serif" )
	  .style( "text-anchor", "middle" )
      .text(function(d) { var s = d.r * 2.4 / fs; if ( d.data.name.length < s ) return d.data.name; else return ""; });
    
  if ( dglabels == "yes" )
  {
		node.filter(function(d) { return d.children && d.data.name != "root"; }).append("text")
			.attr("dy", function( d ) { return ( d.r + fs * 2 ); } )
			.attr("dx", fs )
			.style("font", fs * 1.5 + "px sans-serif" )
			.style("font-weigth", "bolder" )
			.style( "text-anchor", "middle" )
			.text(function(d) { return d.data.name; });
  }


}


draw( "#normal" );

function click(d) {
  console.log( "bb " );
  if (d.title) {
	console.log( " mm " + d.title );
	d._title = d.title;
    d.title = null;
  } else {
	console.log( " nn " + d._title );
    d.title = d._title;
    d._title = null;
  }
  update( d );
}

	d3.select("#generate")
		.on("click", writeDownloadLink);
	d3.select("#gpng")
		.on("click", writeDownloadPng);

function writeDownloadPng(){
		draw( "#highres" );
		var element = document.getElementById('highres');
		element.style.background="white";
		html2canvas(element, {
			onrendered: function(canvas) {
				canvas.toBlob(function(blob) {
					saveAs(blob, outfn);
				}, "image/png");
		} } );
		//draw("#normal" );
	}

	function writeDownloadLink(){
		try {
			var isFileSaverSupported = !!new Blob();
		} catch (e) {
			alert("blob not supported");
		}

		var html = '<svg width="960" height="500" title="test2" version="1.1" xmlns="http://www.w3.org/2000/svg">'
					+ d3.select("#mainsvg").node().innerHTML + '</svg>';

		var blob = new Blob([html], {type: "image/svg+xml"});
		saveAs(blob, "bubblechart.svg");
	};

	
 </script>
 </body>
 </html>
"""