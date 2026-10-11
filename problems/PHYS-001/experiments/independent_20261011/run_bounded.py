"""Run one child on one CPU; preserve binary stdout/stderr and enforce frozen-input hashes."""
import argparse
import ctypes
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--round', required=True, type=Path)
    p.add_argument('--tag', required=True)
    p.add_argument('--replay', action='store_true', help='Explicit reproduction of a finished frozen record, no new novelty claim')
    args = p.parse_args()
    if not args.tag or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789-_' for c in args.tag):
        raise RuntimeError('unsafe output tag')
    record = json.loads(args.round.read_text(encoding='utf-8'))
    root = args.round.resolve().parents[2]
    if not record.get('admitted_at') or record['state'] not in ('ADMITTED','ROUND_FINISHED'):
        raise RuntimeError('missing valid local admission')
    if record['state']=='ROUND_FINISHED' and not args.replay:
        raise RuntimeError('finished round requires explicit --replay')
    for relative, digest in record['execution_inputs_sha256'].items():
        path = root / relative
        if not path.is_relative_to(root) or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise RuntimeError('FAIL_CLOSED: frozen execution input changed: '+relative)
    base = root / 'problems' / record['problem_id'] / 'results' / 'independent_20261011'
    attempt = base / args.tag
    attempt.mkdir(parents=True, exist_ok=False)
    output = attempt / 'data'
    now = dt.datetime.now(dt.timezone.utc)
    elapsed = (now-dt.datetime.fromisoformat(record['started_at'])).total_seconds()
    remaining = 1800 if args.replay else min(1800, record['budget']['wall_minutes']*60-elapsed)
    if remaining<=0:
        (attempt/'launch-error.txt').write_text('Budget exhausted before launch; NOT_RUN.\n',encoding='utf-8')
        return 124
    affinity = {}
    if os.name=='nt':
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.GetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t), ctypes.POINTER(ctypes.c_size_t)]
        kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        proc_mask, system_mask = ctypes.c_size_t(), ctypes.c_size_t()
        current = kernel.GetCurrentProcess()
        if not kernel.GetProcessAffinityMask(current,ctypes.byref(proc_mask),ctypes.byref(system_mask)):
            raise ctypes.WinError(ctypes.get_last_error())
        mask = proc_mask.value & -proc_mask.value
        if not kernel.SetProcessAffinityMask(current,mask):
            raise ctypes.WinError(ctypes.get_last_error())
        affinity = {'original_mask':proc_mask.value, 'selected_mask':mask, 'cpu_count':mask.bit_count(), 'child_inherits_mask':True}
    else:
        allowed = os.sched_getaffinity(0)
        selected = min(allowed)
        os.sched_setaffinity(0,{selected})
        affinity = {'selected_cpu':selected,'cpu_count':1,'child_inherits_mask':True}
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1', PYTHONIOENCODING='utf-8',
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    command = [sys.executable,'-X','utf8',str(HERE/'verify.py'),'--round',str(args.round.resolve()),'--output',str(output)]
    launch = {'started_at':now.isoformat(),'command':command,'cwd':str(root),'python':sys.version,
              'executable':sys.executable,'platform':platform.platform(),'affinity':affinity,
              'timeout_seconds':remaining,'usd':0,'replay':args.replay}
    (attempt/'launch.json').write_text(json.dumps(launch,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'launch':launch}),flush=True)
    t0 = time.perf_counter()
    timed_out = False
    with (attempt/'stdout.txt').open('wb') as stdout, (attempt/'stderr.txt').open('wb') as stderr:
        child = subprocess.Popen(command,cwd=root,env=env,stdout=stdout,stderr=stderr)
        try:
            code = child.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()
            code,timed_out = 124,True
    result_path = output/'result.json'
    summary = json.loads(result_path.read_text(encoding='utf-8')) if result_path.exists() else None
    success = code==0 and summary is not None and summary.get('all_criteria_passed') is True and summary.get('exit_code')==0
    finished = {'finished_at':dt.datetime.now(dt.timezone.utc).isoformat(),'wall_seconds':time.perf_counter()-t0,
                'child_exit_code':code,'timed_out':timed_out,'result_present':summary is not None,
                'execution_status':'COMPLETED_WITH_AUDIT_FAILURE' if code==1 and summary else 'COMPLETED' if success else 'FAILED',
                'scientific_all_criteria_passed':success,'result_verdict':summary.get('verdict') if summary else None}
    (attempt/'execution.json').write_text(json.dumps(finished,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(finished),flush=True)
    return code if code else 0 if success else 1

if __name__=='__main__':
    raise SystemExit(main())
