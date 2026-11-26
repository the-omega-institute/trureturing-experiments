"""Repeat the certified cyclic vacuum with a genuinely CP-conserving CKM target."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb

cyclic = run_path(str(Path(__file__).with_name("ch39-shell-gw-cyclic-vacua.py")))
six = cyclic["six"]
m = six["m"]
arc = m["arc"]
cost = m["m"]
orb = cost["m"]["orb"]
four = cost["m"]["four"]
arr = four["float_matrix"]


def target(delta):
    eye = cost["m"]["eye"]
    v = eye(3)
    for i,j,t,ph in ((1,2,".04275",arb(0)),(0,2,".00377",-delta),(0,1,".22704",arb(0))):
        theta = arb(t)
        b = eye(3)
        b[i,i] = b[j,j] = theta.cos()
        b[i,j] = acb(0,ph).exp()*theta.sin()
        b[j,i] = -b[i,j].conjugate()
        v = v*b
    _,du,dd = orb["reference"]()
    U,A,B = [four["exact_frame"](np.array(z).reshape(3,3)) for z in cost["units"]]
    u = U*du*A*cost["up"][-1]
    d = U*v*dd*B*cost["down"][-1]
    return acb_mat([[z] for z in u.entries()+d.conjugate().entries()])


def seed_for(new_target):
    P = np.vstack([arr(cost["F"]),arr(arc["conjugate_output"](six["t"]))])
    mat = arr(m["K"])
    z = np.array([complex(x,y) for x,y in six["seed"]])
    ref = arr(m["tangent"](acb_mat(54,1),six["q"],six["lam"]))
    goal = arr(new_target).flatten()
    pinv = np.linalg.pinv(P)
    z = z+pinv@(goal-P@z)
    projector = np.eye(54)-pinv@P
    for _ in range(5):
        field = z.reshape(6,9)
        blocks = np.einsum("ir,irjs,js->ij",field.conj(),mat.reshape(6,9,6,9),field)
        g = (blocks+ref).sum(axis=1).imag[:5]
        h = np.einsum("ir,irjs->ijs",field.conj(),mat.reshape(6,9,6,9)).reshape(6,54)
        for i in range(6):
            h[i,9*i:9*i+9] -= (mat@z)[9*i:9*i+9].conj()
        N = projector@h[:5].conj().T
        z -= 1j*N@np.linalg.solve((h[:5]@N).real,g)
    return acb_mat([[acb(float(v.real),float(v.imag))] for v in z])


if __name__=="__main__":
    new_target = target(arb(0))
    z = m["exact_fit"](seed_for(new_target),six["t"],six["lam"],target=new_target)
    six["certify"](z,norm_limit="12.111")
    cyclic["certify"](z)
