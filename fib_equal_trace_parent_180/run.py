#!/usr/bin/env python3
"""Owned, finite experiment entry; leaves mathematical outputs, no journals."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parent
TASKS={
    'certificate': [('same_trace.py',),('verify_same_trace.py',)],
    'fixed-images': [('incidence.py','fixed'),('incidence.py','singular')],
    'universal-trial': [('incidence.py','universal','--modulus','1000003')],
}


def main():
    assert sys.version_info>=(3,12)
    parser=argparse.ArgumentParser()
    parser.add_argument('task',choices=TASKS)
    args=parser.parse_args()
    for script,*arguments in TASKS[args.task]:
        start=time.monotonic()
        # subprocess.run owns, terminates and joins this exact child on timeout.
        try:
            result=subprocess.run([sys.executable,'-B',str(ROOT/script),*arguments],
                                  timeout=60,check=False)
        except subprocess.TimeoutExpired:
            print(json.dumps({'status':'timeout','allowance_seconds':60,
                              'wall_seconds':time.monotonic()-start}),flush=True)
            return 124
        if result.returncode:
            return result.returncode
    return 0


if __name__=='__main__':raise SystemExit(main())
