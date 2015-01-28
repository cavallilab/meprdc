#
# a simple python script which reads a GROMACS style topology file (top) and a datafile in the internal RDC file format and outputs 
# data for an input plumed.dat-file
#
# S OLSSON
#

from sys import argv

pline = lambda i,a,b,d,e:"ATOMS%s=%s,%s COUPLING%s=%s TESTSET%s=%s"%(i,a,b,i,d,i,e)

if __name__=='__main__':
   f_rdcs = open(argv[1]).readlines()
   f_topo = open(argv[2]).readlines()
   for u,line in enumerate(f_rdcs):
      ls=line.split()
      #1 2 N 2 H -1.710
      d=ls[-1]
      i=ls[1]
      i_n=ls[2]
      j=ls[3]
      j_n=ls[4]
      an=[str(u+1),"0","0",d]
      for tline in f_topo:
         if tline[0]==";" or len(tline)<50 or tline[0]=="#":
            continue
         tls=tline.split()
         if tls[2]==i and tls[4]==i_n:
            an[1]=tls[0]
            continue
         
         if tls[2]==j and tls[4]==j_n:
            an[2]=tls[0]
            continue

         if an[1]!="0" and an[2]!="0":
            break
      print pline(an[0], an[1], an[2], an[3],"0")

