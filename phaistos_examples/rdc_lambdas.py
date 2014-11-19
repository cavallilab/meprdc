#
# quick and dirty parsing of inf-rdc output, estimation of Lagrange multipliers via the estimate_scales class
#
#
from estimate_scales import *
from sys import argv
from pdb import *
from copy import deepcopy as dc
if __name__=="__main__":
  base_path="FILL THIS" 

  currn=argv[1]
  ensemble_averages=[]
  calculated_data=[]
  magnetic_fields=[]
  for fn in argv[2:]:
    f=open(fn)
    lines=f.readlines()
    f.close()
    startpos=lines[1].find("inf-rdc{")+len("inf-rdc{")

    for line in lines[1:]:
      i=0
      ea2=[]
      cd2=[]
      magf2=[]
      for alignment in line[startpos:].split('}{'):
        ea=[]
        cd=[]
        magf=[]
        calcd,ensa,mf=alignment.split('];[')
        for u in calcd.strip('](),[(').split('),('):
          cd=cd+map(float, u.split(','))

        for u in ensa.strip('](),[(').split('),('):
          ea=ea+map(float, u.split(','))

        try:
          magf=magf+map(lambda x:float(x.strip(')')), mf.strip('](),[(').replace("]}/","").split(','))
        except:
          magf=magf+map(lambda x:float(x.strip(')')), mf.strip('](),[(').replace("]}/","").split(',')[:-1])
        ea2.append(dc(ea))
        cd2.append(dc(cd))
        magf2.append(dc(magf))

      ensemble_averages.append(dc(ea2))
      calculated_data.append(dc(cd2))
      magnetic_fields.append(dc(magf2))
        
  ensemble_averages=[map(array, s) for s in ensemble_averages]
  data=ensemble_averages[0]
  samples=array([map(array, s) for s in calculated_data]).T


  se=[]
  fff=open(base_path+"/output_ens.dat", "a")
  fff.write(argv[1]+'\n')
  for i in range (len(data)):
    tmp=[aa.tolist() for aa in samples[i]]
    se.append(ScaleEstimator([reshape(data[i], (1,len(data[i])))], [array(tmp).T]))
    se[i].small=1e-6 # Delta t
    fff.write("r:"+str( corrcoef(se[i].data[0][0],se[i].means[0])[0,1])+'\n')
  fff.write('\n')
  fff.close()
  for i in range(10):
    for s in se:
      s.mh_step()
  i=0
  while(i<len(se)): 
    try:
      savetxt(base_path+str(i)+"_"+str(int(currn)), se[i].alla[0][0,:]+loadtxt(base_path+str(i)+'_'+str(int(currn)-1) ))
    except:
      set_trace()

    i=i+1

