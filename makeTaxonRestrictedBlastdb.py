#!/usr/share/python3Venv/orange3Venv/bin/python3
# -*- coding: utf8 -*-

# import clauses
import tempfile, io, sys, argparse, multiprocessing, os, os.path, glob, ete3, numpy, pandas, functools, itertools, pathlib, re


# arguments parser
parser = argparse.ArgumentParser()
parser.add_argument(\
    "--public_nr_blast_dir",action="store",dest='publicBlastDbsDirPath',type=pathlib.Path,required=False,\
    default='/opt/data/public_data/ftp.ncbi.nlm.nih.gov/blast/db',\
    help="  Full path to the directory where to look for the public BLAST nr databases.")
parser.add_argument(\
    "--private_nr_blast_dir",action="store",dest='privateBlastDirPath',type=pathlib.Path,required=False,\
    default='/opt/data/private_data/ftp.ncbi.nlm.nih.gov/blast',\
    help="  Full path to the directory where the taxblast dbs are created.")
parser.add_argument(\
    "--taxid",action='store',dest='taxid',type=int,default=None,required=True,help="  taxid to consider for building the restricted BLAST db.")
parser.add_argument(\
    "--cpu",action='store',dest='cpu',type=int,default=multiprocessing.cpu_count()-4,required=False,help="  number of cpu to use for parallelized steps.")
args=vars(parser.parse_args())
for key in args:exec(key+"=args[\'"+key+"\']")

# enable the use of environment modules
exec(open('/usr/share/modules/init/python.py3','r').read())

# test parameters and define global variables
taxid=str(taxid)

# functions definitions
######################################################################################################
def envModuleLoad(modulefile): module('load',modulefile)

######################################################################################################
def envModuleUnload(modulefile): module('unload',modulefile)

######################################################################################################
def checkTaxid(taxid=None):
  if taxid=='1' or int(taxid) in ete3.NCBITaxa().get_descendant_taxa(1,intermediate_nodes=True): return True
  else: return False

######################################################################################################
def makeTaxidsListFp(taxid=None):
  outFp=privateBlastDirPath/'taxids_gis'/taxid/'taxids')
  if not(outFp.is_file() and outFp.stat().st_size>0):
    outFp.open(mode='wt').writelines(sorted([taxid+'\n']+list(map(lambda i:str(i)+'\n',ete3.NCBITaxa().get_descendant_taxa(taxid)))))
  if not(outFp.is_file() and outFp.stat().st_size>0):return None
  else: return outFp

######################################################################################################
def makeTaxidsGisListFpFromDb(taxid=None,blastDbPrefixPath=None):
  outFp=privateBlastDirPath/'taxids_gis'/taxid/(blastDbPrefixPath.split('/')[-1]+'.taxids_gis')
  if not(outFp.is_file() and outFp.stat().st_size>0):
    with outFp.open('wt') as outFh:
      envModuleLoad('ncbi-blast+/2.5.0')
      blastdbCmdLine="blastdbcmd -db "+blastDbPrefixPath+" -entry all -outfmt \"%T %g\" | sort -k1,1 -"
      a=subprocess.Popen(blastdbCmdLine,shell=True,stdout=outFh,stderr=subprocess.PIPE,universal_newlines=True)
      out,err=a.communicate()
      envModuleUnload('ncbi-blast+/2.5.0')
  if not(outFp.is_file() and outFp.stat().st_size>0): return None
  else: return outFp

######################################################################################################
def makeRootTaxidsGisListFps(cpu=None):
  rootBlastDbPrefixPaths=getBlastDbPrefixPaths(taxid='1')
  rootTaxidsGisListFps=list(privateBlastDirPath/'taxids_gis'/'1'/(path.split('/')[-1]+'.taxids_gis') for path in rootBlastDbPrefixPaths)
  if not(all(fp.is_file() and fp.stat().st_size>0 for fp in rootTaxidsGisListFps)):
    rootTaxidsGisListFps=multiprocessing.Pool(processes=cpu).map(makeTaxidsGisListFpFromDb,list((taxid,path) for path in rootBlastDbPrefixPaths))
  if not(all(fp.is_file() and fp.stat().st_size>0 for fp in rootTaxidsGisListFps)): return None
  else: return rootTaxidsGisListFps
  
######################################################################################################
def getBlastDbPrefixPaths(taxid=None):
  envModuleLoad('ncbi-blast+/2.5.0')
  blastdbcmdCmdLine="blastdbcmd -list_outfmt \"%f\" -list "+str(privateBlastDirPath/'db'/taxid)
  a=subprocess.Popen(blastdbcmdCmdLine,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
  out,err=a.communicate()
  envModuleUnload('ncbi-blast+/2.5.0')
  if out=='': return None
  else: return out.strip().split()

######################################################################################################
def makeGisListFp(taxid=None):
  descGisListFp=privateBlastDirPath/'taxids_gis'/taxid/'gis'
  if not(descGisListFp.is_file() and descGisListFp.stat().st_size>0):
    rootTaxidsGisList=pandas.concat(\
      pandas.read_table(filepath_or_buffer=str(rootTaxidsGisListFp),sep=' ',header=None,names=['taxid','gi'],dtype={'taxid':numpy.int32,'gi': numpy.int64}) \
      for rootTaxidsGisListFp in makeRootTaxidsGisListFps(cpu=cpu))
    descTaxidsList=pandas.read_table(filepath_or_buffer=str(makeTaxidsListFp(taxid=taxid)),header=None,names=['taxid'],dtype={'taxid':numpy.int32})
    descGisList=pandas.merge(rootTaxidsGisList,descTaxidsList,on='taxid',how='inner',sort=False,suffixes=('§gi_list','§taxid_filter'))['gi'].sort_values()
    descGisList.to_csv(path=str(descGisListFp),sep=" ",header=False,index=False)
  if not(descGisListFp.is_file() and descGisListFp.stat().st_size>0): return None
  else: return descGisListFp

######################################################################################################
def makeBlastDbs(taxid=None):
  blastDbAliasFp=privateBlastDirPath/'db'/(taxid+'.nr.pal')
  if not(blastDbAliasFp.is_file() and blastDbAliasFp.stat().st_size>0):
    envModuleLoad('ncbi-blast+/2.5.0')
    blastdbcpCmdLine=\
      'blastdbcp -db '+str(privateBlastDirPath/'db'/'1.nr')+' -gilist '+str(makeGisListFp(taxid=taxid))+\
      ' -out '+str(privateBlastDirPath/'db'/taxid/'nr')+' -membership_bits -title nr'
    a=subprocess.Popen(blastdbcpCmdLine,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    out,err=a.communicate()
    dbList=getBlastDbPrefixPaths(taxid=taxid)
    blastdb_aliastoolCmdLine='blastdb_aliastool -dbtype prot -dblist '+dbList+' -title '+taxid+'.nr -out '+str(privateBlastDirPath/'db'/(taxid+'.nr'))
    a=subprocess.Popen(blastdb_aliastoolCmdLine,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    out,err=a.communicate()
    envModuleUnload('ncbi-blast+/2.5.0')
  if not(blastDbAliasFp.is_file() and blastDbAliasFp.stat().st_size>0): return None
  else: return blastDbAliasFp




dpDict={'root':privateBlastDbsDirPath}
dpDict=makeTaxidDpDict(dpDict=dpDict,privateBlastDbsDirPath=privateBlastDbsDirPath,taxid='1',erase=False)
if checkTaxid(taxid=taxid):
  dpDict=makeTaxidDpDict(dpDict=dpDict,privateBlastDbsDirPath=privateBlastDbsDirPath,taxid=taxid,erase=False)
else: sys.exit('this taxid does nor exist')
taxidsListFps=dict() ; gisListFps=dict() ; blastDbs=dict()
taxidsListFps[taxid]=makeTaxidsListFp(taxid=taxid,dpDict=dpDict)
gisListFps[taxid]=makeGisListFp(taxid=taxid,dpDict=dpDict)
blastDbs[taxid]=makeBlastDbs(taxid=taxid,dpDict=dpDict)

