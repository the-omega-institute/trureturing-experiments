# SPDX-License-Identifier: Apache-2.0
import datetime, json, os, pathlib, signal, subprocess, sys, time
root=pathlib.Path(__file__).resolve().parent
start=time.monotonic(); utc=datetime.datetime.now(datetime.UTC).isoformat()
env=dict(os.environ); env['PATH']=str(pathlib.Path(sys.executable).parent)+':'+env.get('PATH','')+':/usr/sbin'; env['PYTHONDONTWRITEBYTECODE']='1'
pid=os.fork()
if pid==0:
    os.execve(sys.executable,[sys.executable,'-B',str(root/'experiment.py')],env)
peak=0; limit=False; next_sample=start
while True:
    done,status,usage=os.wait4(pid,os.WNOHANG)
    if done: break
    now=time.monotonic()
    if now-start>=90 and not limit:
        os.kill(pid,signal.SIGKILL); limit=True
    if now>=next_sample:
        r=subprocess.run(['/bin/ps','-o','rss=','-p',str(pid)],capture_output=True,text=True)
        if r.stdout.strip().isdigit(): peak=max(peak,int(r.stdout.strip())*1024)
        next_sample=now+0.5
    time.sleep(0.05)
result={'computation_started_utc':utc,'computation_ended_utc':datetime.datetime.now(datetime.UTC).isoformat(),
        'external_hard_limit_seconds':90,'external_limit_triggered':limit,'elapsed_seconds':time.monotonic()-start,
        'exit_code':os.waitstatus_to_exitcode(status),'peak_process_rss_bytes':usage.ru_maxrss,
        'sampled_peak_rss_bytes':peak,'child_user_seconds':usage.ru_utime,'child_system_seconds':usage.ru_stime,
        'own_solver_process_joined':True,'rss_unit':'macOS bytes; ps KiB converted to bytes'}
(root/'resources.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
