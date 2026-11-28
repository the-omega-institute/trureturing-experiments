"""Deterministic complete artifact packing, with bounded Git blob sizes.

Each artifact is either a gzip of its entire logical bytes or those bytes
unchanged, split into <=32 MiB parts. The manifest hashes both levels.
"""
import argparse,gzip,hashlib,json,shutil,tempfile
from pathlib import Path
CHUNK=32*1024*1024

def digest(path):
    h=hashlib.sha256();n=0
    with Path(path).open('rb') as f:
        while block:=f.read(1024*1024):h.update(block);n+=len(block)
    return {'bytes':n,'sha256':h.hexdigest()}

def pack(directory,names):
    directory=Path(directory);mp=directory/'artifacts.json'
    manifest=json.loads(mp.read_text()) if mp.exists() else {'schema':1,'artifacts':{}}
    for name in names:
        source=directory/name;logical=digest(source);parts=[]
        with tempfile.TemporaryDirectory(prefix='pack-',dir=directory) as td:
            packed=Path(td)/'compressed'
            if name.endswith('.gz'):shutil.copyfile(source,packed);encoding='identity'
            else:
                encoding='gzip'
                with source.open('rb') as f,packed.open('wb') as raw:
                    with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0,compresslevel=6) as g:shutil.copyfileobj(f,g,1024*1024)
            compressed=digest(packed)
            with packed.open('rb') as f:
                while block:=f.read(CHUNK):
                    pp=directory/(name+'.part'+str(len(parts)).zfill(3));pp.write_bytes(block);parts.append({'path':pp.name,**digest(pp)})
        manifest['artifacts'][name]={'encoding':encoding,'logical':logical,'packed':compressed,'parts':parts}
        mp.write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
        source.unlink()
    return manifest

def unpack(directory,name,destination):
    directory=Path(directory);m=json.loads((directory/'artifacts.json').read_text())['artifacts'][name]
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='unpack-',dir=destination.parent) as td:
        pp=Path(td)/'packed'
        with pp.open('wb') as out:
            for part in m['parts']:
                p=directory/part['path'];assert digest(p)=={k:part[k] for k in ('bytes','sha256')};
                with p.open('rb') as f:shutil.copyfileobj(f,out)
        assert digest(pp)==m['packed']
        if m['encoding']=='gzip':
            with gzip.open(pp,'rb') as f,destination.open('wb') as out:shutil.copyfileobj(f,out)
        else:shutil.copyfile(pp,destination)
    assert digest(destination)==m['logical'];return m['logical']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--directory',required=True,type=Path);ap.add_argument('--pack',nargs='+');ap.add_argument('--unpack');ap.add_argument('--destination',type=Path);a=ap.parse_args()
    print(json.dumps(pack(a.directory,a.pack) if a.pack else unpack(a.directory,a.unpack,a.destination),sort_keys=True))
