"""Independent direct DFT against reviewed FFT/energy definitions, small N only."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
TOTAL_CHECKS = ('components_and_mutations_pass','case_budget_pass','historical_FAIL_preserved')


def complete(checks, expected):
    return bool(expected) and set(checks) == set(expected) and all(checks[k] is True for k in expected)


def direct_dft(u):
    n = len(u)
    return [complex(math.fsum(v * math.cos(2*math.pi*j*k/n) for j,v in enumerate(u)),
                    -math.fsum(v * math.sin(2*math.pi*j*k/n) for j,v in enumerate(u)))
            for k in range(n)]


def fixture(n, dc, nyquist, modes):
    u = [dc + nyquist*((-1)**j) + math.fsum(a*math.cos(2*math.pi*k*j/n+phase)
                                           for k,a,phase in modes) for j in range(n)]
    spectrum = [0j]*n
    spectrum[0] = complex(n*dc,0)
    spectrum[n//2] = complex(n*nyquist,0)
    for k,a,phase in modes:
        coefficient = n*a/2*complex(math.cos(phase),math.sin(phase))
        spectrum[k] += coefficient
        spectrum[-k] += coefficient.conjugate()
    energy = dc*dc + nyquist*nyquist + math.fsum(a*a/2 for k,a,phase in modes)
    return u, spectrum, energy


def coefficient_error(a,b):
    return max((abs(x-y) for x,y in zip(a,b)), default=math.inf) if len(a)==len(b) else math.inf


def reviewed_functions(name,n):
    path = HERE / 'reference' / (name + '.py')
    wanted = {'fft','energy_budget'}
    nodes = [node for node in ast.parse(path.read_text(encoding='utf-8')).body
             if isinstance(node,ast.FunctionDef) and node.name in wanted]
    if {node.name for node in nodes} != wanted:
        raise RuntimeError('FAIL_CLOSED: reviewed function missing')
    namespace = {'math':math,'N':n}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),namespace)
    return namespace['fft'],namespace['energy_budget']


def pairs(values):
    return [[v.real,v.imag] for v in values]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--round',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    record=json.loads(args.round.read_text(encoding='utf-8'))
    proto=json.loads((HERE/'PROTOCOL.json').read_text(encoding='utf-8'))
    protocol_hash=hashlib.sha256((HERE/'PROTOCOL.json').read_bytes()).hexdigest()
    if record.get('state') not in ('ADMITTED','ROUND_FINISHED') or not record.get('admitted_at') or record.get('protocol_sha256')!=protocol_hash:
        raise RuntimeError('FAIL_CLOSED: missing admission or changed frozen protocol')
    manifest=json.loads((HERE/'SOURCE_MANIFEST.json').read_text(encoding='utf-8'))
    for entry in manifest:
        if hashlib.sha256((HERE/entry['snapshot']).read_bytes()).hexdigest()!=entry['sha256']:
            raise RuntimeError('FAIL_CLOSED: source hash mismatch')
    started,cpu_started=time.perf_counter(),time.process_time()
    args.output.mkdir(parents=True,exist_ok=False)
    coeff_tol=proto['tolerances']['max_absolute_complex_coefficient_error']
    energy_tol=proto['tolerances']['absolute_energy_and_parseval_error']
    relative_tol=proto['tolerances']['original_relative_energy_tolerance']
    modes=[(mode['k'],mode['amplitude'],mode['phase_radians']) for mode in proto['signal']['mixed']['modes']]
    specs=[('zero',0,0,[]),('dc_0.25',0.25,0,[]),('cos_k3_A1.25_phase0',0,0,[(3,1.25,0)]),
           ('sin_k5_A0.75',0,0,[(5,0.75,-math.pi/2)]),('nyquist_A0.5',0,0.5,[]),
           ('mixed',proto['signal']['mixed']['dc'],proto['signal']['mixed']['nyquist'],modes)]
    cases=[]
    expected_names=['N'+str(n)+'_'+name for n in proto['signal']['N']
                    for name in proto['signal']['fixtures']+['energy_clean_and_1.5x']
                    +['reject_'+name for name in proto['mutations']]]
    expected_names += ['summary_failed_check','summary_missing_check','summary_empty_checks']
    samples=[]
    def case(name,passed,**details):
        cases.append({'name':name,'passed':bool(passed),**details})
        print(json.dumps(cases[-1],sort_keys=True),flush=True)
    maxima={'coefficient_error':0.0,'parseval_error':0.0,'analytic_energy_error':0.0}
    for n in proto['signal']['N']:
        fft1,energy1=reviewed_functions('r1',n)
        fft2,energy2=reviewed_functions('r2',n)
        mixed_u=mixed_dft=None
        for name,dc,nyquist,mode_spec in specs:
            u,analytic_spectrum,analytic_energy=fixture(n,dc,nyquist,mode_spec)
            if name=='sin_k5_A0.75':
                # Generate the requested sine directly, not by reusing the analytic phase expression.
                u=[0.75*math.sin(2*math.pi*5*j/n) for j in range(n)]
                analytic_spectrum=[0j]*n
                analytic_spectrum[5]=complex(0,-n*0.75/2)
                analytic_spectrum[-5]=complex(0,n*0.75/2)
            dft=direct_dft(u)
            fft,fft_r2=fft1(u),fft2(u)
            mean_square=math.fsum(v*v for v in u)/n
            parseval_dft=math.fsum(abs(c)**2 for c in dft)/(n*n)
            parseval_fft=math.fsum(abs(c)**2 for c in fft)/(n*n)
            metrics={'dft_analytic_coefficient_error':coefficient_error(dft,analytic_spectrum),
                     'fft_dft_coefficient_error':coefficient_error(fft,dft),
                     'fft_analytic_coefficient_error':coefficient_error(fft,analytic_spectrum),
                     'r2_fft_dft_coefficient_error':coefficient_error(fft_r2,dft),
                     'dft_parseval_error':abs(mean_square-parseval_dft),
                     'fft_parseval_error':abs(mean_square-parseval_fft),
                     'analytic_energy_error':abs(mean_square-analytic_energy)}
            maxima['coefficient_error']=max(maxima['coefficient_error'],*(v for k,v in metrics.items() if 'coefficient' in k))
            maxima['parseval_error']=max(maxima['parseval_error'],metrics['dft_parseval_error'],metrics['fft_parseval_error'])
            maxima['analytic_energy_error']=max(maxima['analytic_energy_error'],metrics['analytic_energy_error'])
            valid=all(math.isfinite(v) and v <= (coeff_tol if 'coefficient' in k else energy_tol) for k,v in metrics.items())
            case('N'+str(n)+'_'+name,valid,metrics=metrics)
            samples.append({'N':n,'fixture':name,'input':u,'analytic_spectrum':pairs(analytic_spectrum),
                            'dft':pairs(dft),'fft_r1':pairs(fft),'fft_r2':pairs(fft_r2),
                            'mean_square':mean_square,'analytic_energy':analytic_energy,
                            'parseval_dft':parseval_dft,'parseval_fft':parseval_fft})
            if name=='mixed':
                mixed_u,mixed_dft=u,dft
        energy_u,_,energy_analytic=fixture(n,0,0,modes)
        amps=[a for k,a,phase in modes]
        inflated=[1.5*v for v in energy_u]
        clean_residual,dirty_residual=energy1(energy_u,amps),energy1(inflated,amps)
        r2_clean,r2_dirty=energy2(energy_u,amps),energy2(inflated,amps)
        energy_metrics={'reference_residual':clean_residual,'inflated_residual':dirty_residual,
                        'r2_component_reference_residual':r2_clean,'r2_component_inflated_residual':r2_dirty}
        name='N'+str(n)+'_energy_clean_and_1.5x'
        case(name,clean_residual<=relative_tol and dirty_residual>relative_tol
             and r2_clean<=relative_tol and r2_dirty>relative_tol
             and abs(dirty_residual-1.25)<=energy_tol,metrics=energy_metrics)
        sign_error=coefficient_error([c.conjugate() for c in mixed_dft],mixed_dft)
        normalized_error=coefficient_error([c/n for c in mixed_dft],mixed_dft)
        spectrum_sum=math.fsum(abs(c)**2 for c in mixed_dft)
        mixed_mean=math.fsum(v*v for v in mixed_u)/n
        wrong_parseval_N=abs(mixed_mean-spectrum_sum/n)
        wrong_parseval_N3=abs(mixed_mean-spectrum_sum/(n**3))
        poisoned_reference=energy1(inflated,[1.5*a for a in amps])
        wrong_single_reference=energy1(energy_u,[1.0])
        mutations=[('positive_exponent_conjugation',sign_error>coeff_tol,{'coefficient_error':sign_error,
                      'parseval_error_still_small':abs(mixed_mean-spectrum_sum/(n*n))}),
                   ('forward_divide_by_N',normalized_error>coeff_tol,{'coefficient_error':normalized_error}),
                   ('Parseval_divide_by_N_instead_of_N_squared',wrong_parseval_N>energy_tol,{'energy_error':wrong_parseval_N}),
                   ('Parseval_divide_by_N_cubed',wrong_parseval_N3>energy_tol,{'energy_error':wrong_parseval_N3}),
                   ('inflated_reference_amplitudes_1.5',poisoned_reference<=relative_tol,
                    {'inflated_input_poisoned_reference_residual':poisoned_reference,
                     'mutation_guard_incorrectly_rejects_inflated':poisoned_reference>relative_tol}),
                   ('wrong_single_amplitude_reference_[1.0]',wrong_single_reference>relative_tol,
                    {'clean_wrong_reference_residual':wrong_single_reference,
                     'original_guard_form_is_discriminating':False})]
        for mutation,rejected,details in mutations:
            name='N'+str(n)+'_reject_'+mutation
            case(name,rejected,details=details)
        samples.append({'N':n,'fixture':'energy_only','input':energy_u,'amplitudes':amps,
                        'inflated_input':inflated,'analytic_energy':energy_analytic,**energy_metrics})
    for name,flags,expected in [('summary_failed_check',{'a':True,'b':False},['a','b']),
                                ('summary_missing_check',{'a':True},['a','b']),('summary_empty_checks',{},[])]:
        case(name,complete(flags,expected) is False)
    inherited=json.loads((HERE/'reference/inherited-r1.json').read_text(encoding='utf-8'))
    guard=json.loads((HERE/'reference/inherited-r2-guard.json').read_text(encoding='utf-8'))
    negative=next(c for c in inherited['cases'] if c['beta_name']=='steep_negative')
    preserved={'r1_E3_mean':negative['sf_mean_slope'],'r1_E3_threshold':inherited['config']['tolerances']['neg_min'],
               'r1_E3_check':inherited['checks']['E3_negative_control_steep'],'r1_verdict':inherited['verdict'],
               'r2_amplitude_guard_historical':guard,'classification':'INHERITED_NOT_RERUN'}
    inherited_pass=preserved['r1_E3_check'] is False and preserved['r1_E3_mean']<preserved['r1_E3_threshold']
    inherited_pass=inherited_pass and preserved['r1_verdict']=='BASELINE_FAIL_DO_NOT_PROCEED_TO_NOVELTY'
    flags={c['name']:c['passed'] for c in cases}
    checks={'components_and_mutations_pass':complete(flags,expected_names),
            'case_budget_pass':len(cases)<=100,'historical_FAIL_preserved':inherited_pass}
    passed=complete(checks,TOTAL_CHECKS)
    result={'protocol_sha256':protocol_hash,'checks':checks,'all_criteria_passed':passed,
            'verdict':'SYNTHETIC_COMPONENT_AUDIT_PASS' if passed else 'SYNTHETIC_COMPONENT_AUDIT_FAIL',
            'exit_code':0 if passed else 1,'max_errors':maxima,'named_cases':cases,'historical':preserved,
            'wall_seconds':time.perf_counter()-started,'cpu_seconds':time.process_time()-cpu_started,
            'not_run':proto['not_run'],'physical_discovery_claim':False,
            'independence':'direct O(N^2) DFT and analytic coefficients vs original radix-2 FFT; source-read, not blinded'}
    (args.output/'samples.json').write_text(json.dumps(samples,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    (args.output/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps({'summary':checks,'verdict':result['verdict'],'max_errors':maxima,'exit_code':result['exit_code']},sort_keys=True),flush=True)
    return result['exit_code']

if __name__=='__main__':
    raise SystemExit(main())
