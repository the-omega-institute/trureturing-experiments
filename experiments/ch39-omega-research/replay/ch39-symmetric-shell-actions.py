"""Flux and face Maxwell action on all coordinate-permutation invariant convex shell boundaries."""
from itertools import combinations, permutations
from pathlib import Path
import runpy

from flint import arb

geometry = runpy.run_path(str(Path(__file__).with_name('ch39-shell-chiral-flux.py')))
sites = geometry['pts']
edges = geometry['edges']
negative = {pair for pair, (distance, _) in edges.items() if distance == 200}

def orbits(symmetries):
    symmetries = tuple(symmetries)
    return sorted({tuple(sorted({sites.index(tuple(site[i] for i in permutation))
                                for permutation in symmetries})) for site in sites})

def area(face):
    normal = [0, 0, 0]
    for i, j in zip(face, face[1:] + face[:1]):
        term = geometry['cross'](sites[i], sites[j])
        normal = [a+b for a,b in zip(normal,term)]
    return arb(sum(x*x for x in normal)).sqrt()/2

def candidates(symmetries):
    blocks = orbits(symmetries)
    result = []
    for count in range(1, len(blocks)+1):
        for chosen in combinations(blocks,count):
            vertices = tuple(sorted(x for block in chosen for x in block))
            faces = geometry['boundary'](vertices)
            if not faces or any(tuple(sorted((i,j))) not in edges for face in faces
                                for i,j in zip(face,face[1:]+face[:1])):
                continue
            degree, margin = geometry['degree'](faces,negative)
            areas = [area(face) for face in faces]
            fluxes = [geometry['lifted'](face,negative) for face in faces]
            action = sum(phi**2/a for phi,a in zip(fluxes,areas))
            result.append((vertices,degree,sum(areas),action))
    return result

if __name__=='__main__':
    full = candidates(tuple(permutations(range(3))))
    cyclic = candidates(((0,1,2),(1,2,0),(2,0,1)))
    small = geometry['inner']
    zero = (11,12,14,18)
    symmetric = {row[0]:row for row in full}
    first, second = symmetric[small], symmetric[zero]
    assert len(full)==10 and len(cyclic)==22 and first[1]==1 and second[1]==0
    critical = (first[3]-second[3])/(second[2]-first[2])
    assert arb('.00616333')<critical<arb('.00616334')
    rivals = ((9,14,15,16),(10,13,14,17))
    # The cyclic competitors have the same areas by congruence. Check their energies
    # agree exactly by the edge phases of their three triangular sides and base.
    def face_signature(vertices):
        signature = []
        for face in geometry['boundary'](vertices):
            chords = sorted(edges[tuple(sorted((i,j)))][0] for i,j in zip(face,face[1:]+face[:1]))
            a,b,c = (sites[i] for i in face)
            triple = abs(geometry['dot'](a,geometry['cross'](b,c)))
            signature.append((tuple(chords),triple))
        return sorted(signature)
    for row in cyclic:
        if row[0] in (small,zero):
            continue
        if row[0] in rivals:
            assert face_signature(row[0])==face_signature(zero) and row[1]==0
        else:
            assert row[2]>second[2] and row[3]>second[3]
    wall = symmetric[geometry['wall']]
    assert wall[2]>first[2] and wall[3]>first[3]
    print('symmetric and cyclic counts',len(full),len(cyclic))
    print('tension at transition',critical.str(15))
    print('degree-three mirror versus inner tetrahedron: Maxwell and area differences',
          (wall[3]-first[3]).str(12),(wall[2]-first[2]).str(12))
