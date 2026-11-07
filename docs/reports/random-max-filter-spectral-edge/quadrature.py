# SPDX-License-Identifier: Apache-2.0
"""Deterministic beta/Jacobi quadrature; numerical screening, not exact evidence."""
import json
from functools import lru_cache
import numpy as np
from scipy.special import roots_jacobi, ndtri
from scipy.linalg import eigvalsh

@lru_cache(maxsize=None)
def beta_nodes(a, b, nodes):
    x, w = roots_jacobi(nodes, b - 1, a - 1)
    return (x + 1) / 2, w / w.sum()

def matrix(n, nodes, half):
    q = (lambda u: ndtri((1 + u) / 2)) if half else ndtri
    out = np.empty((n, n))
    for j in range(1, n + 1):
        v, wv = beta_nodes(j, n + 1 - j, nodes)
        qv = q(v)
        out[j-1, j-1] = np.dot(wv, qv*qv)
        for i in range(1, j):
            t, wt = beta_nodes(i, j-i, nodes)
            value = np.dot(wv*qv, q(v[:,None]*t[None,:]) @ wt)
            out[i-1, j-1] = out[j-1, i-1] = value
    return out

def main():
    results=[]
    print('group,d,nodes,d_squared_lambda_min,target,absolute_node_doubling_change,relative_trace_error')
    for n in [4, 8, 16, 32, 64, 96]:
        for half in [False, True]:
            values=[]
            for nodes in [220, 440, 880, 1760, 3520]:
                m=matrix(n,nodes,half)
                edge=n*n*eigvalsh(m,subset_by_index=[0,0])[0]
                values.append((edge,abs(np.trace(m)/n-1)))
                if len(values)>1 and abs(values[-1][0]-values[-2][0]) < 1e-5:
                    break
            change=abs(values[-2][0]-values[-1][0])
            target=np.pi/(8 if half else 2)
            row={'group':'B' if half else 'S','d':n,'nodes':nodes,'d_squared_lambda_min':values[-1][0], 'target':target,'node_doubling_change':change,'relative_trace_error':values[-1][1]}
            results.append(row)
            print(','.join(str(row[k]) for k in row),flush=True)
            assert change < 1e-5, (n,half,'quadrature not converged',change)
            assert values[-1][1] < 4e-5, (n,half,'trace not converged')
            assert target*n*n/((n+1)*(n+2))-2e-5 <= values[-1][0] <= target+0.03
    for group in ['S','B']:
        rows=[r for r in results if r['group']==group]
        assert abs(rows[-1]['d_squared_lambda_min']-rows[-1]['target']) < 0.005
        assert abs(rows[-1]['d_squared_lambda_min']-rows[-1]['target']) < abs(rows[1]['d_squared_lambda_min']-rows[1]['target'])
    print('QUADRATURE_CHECKS_OK')
    with open('quadrature-results.json','w') as f:json.dump(results,f,indent=2);f.write('\n')

if __name__=='__main__': main()
