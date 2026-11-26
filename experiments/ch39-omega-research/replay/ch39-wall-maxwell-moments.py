"""Numerical P1 finite elements for degree-three zero-mode moments on the flat shell."""
from pathlib import Path
import runpy
import numpy as np

work = Path(__file__).resolve().parent
geometry = runpy.run_path(str(work / "ch39-symmetric-shell-actions.py"))
uniform = runpy.run_path(str(work / "ch39-wall-intrinsic-coordinate.py"))
alpha, beta, gamma = (float(uniform[k]) for k in ("alpha", "beta", "gamma"))
pole = np.array([0., 0.])
mid = np.array([np.sqrt(18), 0.])
end = mid + [2, 0]
inner = mid + [0, 5*np.sqrt(2)]
outer = end + [0, 7*np.sqrt(2)]
south = end + [7*np.sqrt(2/3), 0]
corners = np.array([pole, inner, outer, south])
centroid = mid + [19/18, 0]
vertices = np.array([pole, mid, end, south, inner, outer])
pieces = [(0,1,4), (1,2,4), (4,2,5), (2,3,5)]
faces = geometry["geometry"]
fields = {}
for face in faces["big"]:
    tag = "middle" if len(face) == 4 else "top" if 14 in face else "base"
    fields[tag] = float(faces["lifted"](face, geometry["negative"])/geometry["area"](face))
quadrature_nodes, quadrature_weights = np.polynomial.legendre.leggauss(6)
quadrature_nodes = (quadrature_nodes + 1)/2
quadrature_weights /= 2

# The two columns give D0 = log(H/lambda) and U = log(|z|/r).
singular = np.array([[1/(3*alpha)-1, 1/(3*alpha)], [1/beta-1, 0],
                     [1/gamma-1, 0], [4, -1]])


def mesh(n):
    ids, positions, triangles, magnetic = {}, [], [], []
    for (a,b,c), layer in zip(pieces, ("top","middle","middle","base")):
        local = {}
        for i in range(n+1):
            for j in range(n-i+1):
                key = tuple(sorted((k,s) for k,s in ((a,n-i-j),(b,i),(c,j)) if s))
                if key not in ids:
                    ids[key] = len(positions)
                    positions.append(sum(vertices[k]*s for k,s in key)/n)
                local[i,j] = ids[key]
        for i in range(n):
            for j in range(n-i):
                triangles.append((local[i,j],local[i+1,j],local[i,j+1]))
                magnetic.append(fields[layer])
                if i+j+1 < n:
                    triangles.append((local[i+1,j],local[i+1,j+1],local[i,j+1]))
                    magnetic.append(fields[layer])
    return np.asarray(positions), np.asarray(triangles), np.asarray(magnetic)


def solve_neumann(weights, rows, cols, rhs, anchor):
    diagonal = np.bincount(rows[rows == cols], weights=weights[rows == cols])
    # Pin one nodal value. A constant can be restored after the solve.
    def apply(v):
        result = np.bincount(rows, weights=weights*v[cols], minlength=len(v))
        result[anchor] = 0
        return result
    vector = np.zeros_like(rhs)
    residual = rhs.copy()
    residual[anchor] = 0
    direction = residual/diagonal
    magnitude = residual @ direction
    tolerance = 1e-11*np.linalg.norm(residual)
    for _ in range(12*len(rhs)):
        product = apply(direction)
        step = magnitude/(direction @ product)
        vector += step*direction
        residual -= step*product
        if np.linalg.norm(residual) < tolerance:
            return vector
        z = residual/diagonal
        updated = residual @ z
        direction = z + updated/magnitude*direction
        magnitude = updated
    raise ArithmeticError("finite element linear solve did not converge")


def compute(n):
    points, triangles, magnetic = mesh(n)
    coords = points[triangles]
    a,b,c = coords.transpose(1,0,2)
    twice_area = (b[:,0]-a[:,0])*(c[:,1]-a[:,1])-(b[:,1]-a[:,1])*(c[:,0]-a[:,0])
    assert np.all(twice_area>0)
    edges = np.roll(coords,-1,axis=1)-np.roll(coords,-2,axis=1)
    stiffness = np.einsum("nik,njk->nij", edges, edges)/(2*twice_area[:,None,None])
    rows, cols = np.broadcast_arrays(triangles[:,:,None],triangles[:,None,:])
    rhs = np.zeros((len(points),2))
    for j in range(3):
        np.add.at(rhs[:,0],triangles[:,j],magnetic*twice_area/3)
    boundary = {}
    for tri in triangles:
        for i,j in zip(tri,np.roll(tri,-1)):
            if (j,i) in boundary:
                del boundary[j,i]
            else:
                boundary[i,j] = True
    # Boundary normal derivative of the regular part cancels the explicit log poles.
    for i,j in boundary:
        side = points[j]-points[i]
        normal = np.array([side[1],-side[0]])
        top = (points[i]-corners) @ normal
        top[np.abs(top)<1e-11] = 0
        for t,w in zip(quadrature_nodes,quadrature_weights):
            squared = np.sum((points[i]+t*side-corners)**2,axis=1)
            flux = (top/squared) @ singular
            rhs[i] -= w*(1-t)*flux
            rhs[j] -= w*t*flux
    assert abs(rhs.sum(axis=0)).max()<2e-8
    anchor = int(np.argmin(np.linalg.norm(points-centroid,axis=1)))
    reg = np.column_stack([solve_neumann(stiffness.ravel(), rows.ravel(),cols.ravel(),
                                         rhs[:,i],anchor) for i in range(2)])
    boundary_nodes = np.flatnonzero((np.abs(points[:,1])<1e-10)
                                    & (points[:,0]>=mid[0]-1e-10)
                                    & (points[:,0]<=end[0]+1e-10))
    order = boundary_nodes[np.argsort(points[boundary_nodes,0])]
    reg[:,1] -= (np.interp(centroid[0],points[order,0],reg[order,1])
                 + np.log(np.linalg.norm(centroid-corners,axis=1)) @ singular[:,1])
    moments = np.zeros(3)
    for s,ws in zip(quadrature_nodes,quadrature_weights):
        for t,wt in zip(quadrature_nodes,quadrature_weights):
            bary = np.array([1-s,s*(1-t),s*t])
            samples = np.einsum("ntk,t->nk",coords,bary)
            logs = np.log(np.linalg.norm(samples[:,None,:]-corners,axis=2)) @ singular
            values = logs + np.einsum("ntk,t->nk",reg[triangles],bary)
            moments += (twice_area*ws*wt*s) @ np.exp(values[:,0:1]+2*values[:,1:2]*np.arange(3))
    moments /= moments.sum()
    dual = 1/moments
    dual /= dual.sum()
    omega = np.exp(2j*np.pi/3)
    return moments, abs(dual @ np.array([1,omega,omega**2]))


if __name__ == "__main__":
    for n in (12,24,48):
        moments, coherence = compute(n)
        print("subdivisions, normalized moments, scalar coherence",n,moments,coherence)
