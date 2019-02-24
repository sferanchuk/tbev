library(bio3d)
pdb <- read.pdb("pri-znmg-430-480-f.pdb")
modes <- nma(pdb)
print( modes )
dput( modes$modes[,7], file="modes7.r" )
dput( modes$modes[,8], file="modes8.r" )
dput( modes$modes[,9], file="modes9.r" )


