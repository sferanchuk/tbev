library(bio3d)
args = commandArgs(trailingOnly=TRUE)
pdb <- read.pdb( args[1] )
modes <- nma(pdb)
dim = strtoi( args[2] )
#print( modes )
#dput( modes$modes[,args[2] ], file=args[3] )
dput( modes$modes[,dim] )

