# DTAMefficiency

```bash
cmsrel CMSSW_16_0_0
cd CMSSW_16_0_0/src
cmsenv

git clone git@github.com:battibass/DTNtuples.git DTDPGAnalysis/DTNtuples

cd DTDPGAnalysis/DTNtuples/

git clone git@github.com:CarlosGarciaSanchez19/DTAMefficiency.git

cp DTAMefficiency/* test/
cd test/



```



```bash
# $1 is the ntuples directory, $2 the sample (it must match the name of the ntuple file)
# $3 must be a python list with the desired qualities (for example "['All','nothreehits']" or "['All']" double commas and single commas included)
# $4 is either -a, -nt or -p as well as $5 and $6
# if -a is included the eff is computed using the acceptance
# if -nt is included time cut is not applied (t0 <= 16 TDCs for AM and t0 <= 12.5 ns for RECO)
# if -p is included the efficiency plots are created and saved as pdfs

python3 runAMEffs.py -n $4 $5 $6 -d $1 -s $2 -q $3 
```
Example of a quick analysis and creation of efficiency plots:
```bash

python3 runAMEffs.py -a -p -d /afs/cern.ch/user/c/cgarcias/public/ntuples/AMv1 -s pg-c0 -q "['nothreehits']" 

```
