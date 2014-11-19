from sys import argv
base_path = "BASE PATH HERE"
data_filename = "DATA FILENAME HERE IN BASE PATH"
degree_of_alignment = 0.001
steps_per_sweep = 500000
num_of_alignments = 4

val=str(int(argv[1])-1)
template="""
pdb-file = \"%BASE_PATH%/native.pdb\"
init-from-pdb = 1

iterations = %SPS% 
monte-carlo-metropolis-hastings = 1 


move-none = 1
move-crisp-dbn-eh = 1
move-sidechain-uniform = 1

energy-inf-rdc = 1 # [=arg(= )]                 Activate energy-inf-rdc [number of occurrences]
energy-inf-rdc-weight = 1.
energy-inf-rdc-rdcs-filename = \"%BASE_PATH%/%DATA_FILENAME%\"  
energy-inf-rdc-tensor-model = \"random_magnetic_field\" 
energy-inf-rdc-likelihood = log_linear_ensemble:[scale-matrices:\"%SCLMATS%\"]
energy-inf-rdc-degree-of-alignment=%DOA%


observable-inf-rdc = 1 
observable-inf-rdc-register-interval=2500  
observable-inf-rdc-output-target = "%BASE_PATH%/%k/observables_%t_%n.dat" 
observable-pdb = 1
observable-pdb-register-interval=5000  
observable-pdb-output-target = "%BASE_PATH%/%k/samples_%p_%t_%n/sample_%i_%e.pdb"
"""
template=template.replace("%i",val)
template=template.replace("%n",argv[2])
template=template.replace("%k",argv[1])
template=template.replace("%SCLMATS%", '|'.join(["%BASE_PATH%/"+str(ac)+"_%i" for ac in range(num_of_alignments) ]) )
template=template.replace("%BASE_PATH%", base_path)
template=template.replace("%DATA_FILENAME%", data_filename)
template=template.replace("%DOA%", str(degree_of_alignment))
template=template.replace("%SPS%", str(steps_per_sweep))
opfn=argv[3]
f=open(opfn,"w")
f.write(template)
f.close()
