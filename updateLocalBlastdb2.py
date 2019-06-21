#!/opt/pythonVirtualEnvs/3.5.2/bin/python3

# -*- coding: utf8 -*-

# import clauses
import ftplib, subprocess, filecmp, tempfile, fnmatch, functools, itertools, glob, sys, os, tarfile, multiprocessing, argparse, pathlib

# enable the use of environment modules
exec(open('/usr/share/modules/init/python.py','r').read())

# definition of global variables
blastDbArchivesUri='ftp.ncbi.nlm.nih.gov/blast/db'
localBlastDbArchivesDirPath=pathlib.Path('/opt/data/public_data')/blastDbArchivesUri
localBlastDbDirPath=pathlib.Path('/opt/data/private_data')/blastDbArchivesUri
if not localBlastDbArchivesDirPath.is_dir(): localBlastDbArchivesDirPath.mkdir(parents=True)

if not localBlastDbDirPath.is_dir(): localBlastDbDirPath.mkdir(parents=True)

# functions definitions
######################################################################################################
def envModuleLoad(modulefile): module('load',modulefile)

######################################################################################################
def envModuleUnload(modulefile): module('unload',modulefile)

######################################################################################################
def updateMd5s(localBlastDbArchivesDirPath=None):
  if localBlastDbArchivesDirPath == None or not localBlastDbArchivesDirPath.is_dir(): return None
  with ftplib.FTP('ftp.ncbi.nlm.nih.gov') as ftp:
    ftp.login(user='anonymous') ; ftp.cwd('blast/db')
    remoteMd5FilePathes=list(map(lambda fn:'blast/db/'+fn,fnmatch.filter(dict(ftp.mlsd()).keys(),'nr.*.tar.gz.md5')))
  with pathlib.Path(tempfile.TemporaryDirectory().name) as tmpDir:
    envModuleLoad('aspera-ascp/3.5.4')
    ascpCmdLine="ascp-key -d -T -k 1 --mode=recv --user=anonftp --host=ftp.ncbi.nlm.nih.gov "+' '.join(remoteMd5FilePathes)+" "+str(tmpDir)
    a=subprocess.Popen(ascpCmdLine,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    out, err = a.communicate()
    envModuleUnload('aspera-ascp/3.5.4')
    remoteMd5s=dict((fp.name,fp.open('r').read()) for fp in list(tmpDir.glob('nr.*.tar.gz.md5')))
  localMd5s=dict((fp.name,fp.open('r').read()) for fp in list(localBlastDbArchivesDirPath.glob('nr.*.tar.gz.md5')))
  toDropMd5Dict=dict((fn,fc) for fn,fc in localMd5s.items() if (fn,fc) not in remoteMd5s.items())
  toDownloadMd5Dict=dict((fn,fc) for fn,fc in remoteMd5s.items() if (fn,fc) not in localMd5s.items())
  for fp in map(lambda fn:localBlastDbArchivesDirPath/fn,toDropMd5Dict.keys()): fp.unlink()
  for fn,fc in toDownloadMd5Dict.items(): (localBlastDbArchivesDirPath/fn).open('w').write(fc)
  return {
    'toDrop':list('.'.join(fn.split('.')[:-3]) for fn in toDropMd5Dict.keys()),
    'toDownload':list('.'.join(fn.split('.')[:-3]) for fn in toDownloadMd5Dict.keys())}

######################################################################################################
def dropBlastdbs(localBlastDbArchivesDirPath=None,localBlastDbDirPath=None,dbNamesList=None):
  fpList=list(itertools.chain(*(list(localBlastDbDirPath.glob(dbName+'*'))+[localBlastDbArchivesDirPath/(dbName+'.tar.gz')] for dbName in dbNamesList)))
  for fp in fpList: fp.unlink()
  return not any(fp.is_file() for fp in fpList)

######################################################################################################
def downloadBlastdb(localBlastDbArchivesDirPath=None,localBlastDbDirPath=None,dbName=None):
  envModuleLoad('aspera-ascp/3.5.4')
  print(dbName)
  if not (localBlastDbArchivesDirPath/(dbName+'.tar.gz')).is_file():
    ascpCmdLine="ascp-key -d -T -k 1 --mode=recv --user=anonftp --host=ftp.ncbi.nlm.nih.gov blast/db/"+dbName+'.tar.gz '+str(localBlastDbArchivesDirPath)
    a=subprocess.Popen(ascpCmdLine,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    out=a.communicate()
    envModuleUnload('aspera-ascp/3.5.4')
  if len(list(localBlastDbDirPath.glob(dbName+'*')))==0:
    tarfile.open(str(localBlastDbArchivesDirPath/(dbName+'tar.gz')),'r|gz').extractall(path=str(localBlastDbDirPath))
  return \
    (lambda ll:ll if all(fp.is_file and fp.stat().st_size>0 for fp in ll) else None)\
    ([localBlastDbArchivesDirPath/(dbName+'.tar.gz')]+list(localBlastDbDirPath.glob(dbName+'*')))

updateDict=updateMd5s(localBlastDbArchivesDirPath=localBlastDbArchivesDirPath)
print(updateDict)
dropBlastdbs(localBlastDbArchivesDirPath=localBlastDbArchivesDirPath,localBlastDbDirPath=localBlastDbDirPath,dbNamesList=updateDict['toDrop'])

if __name__ == '__main__':
    with multiprocessing.Pool(processes=24) as pool:
        print(pool.starmap(downloadBlastdb,list([localBlastDbArchivesDirPath,localBlastDbDirPath,dbName] for dbName in updateDict['toDownload'])))
