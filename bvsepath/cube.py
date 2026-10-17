from pathlib import Path
import numpy as np

BOHR_TO_ANG = 0.529177210903


class Cube:
    def __init__(self, comment1, comment2, origin_A, shape, grid_A, cell_A, atoms, values):
        self.comment1 = str(comment1)
        self.comment2 = str(comment2)
        self.origin_A = np.asarray(origin_A, float)
        self.shape = np.asarray(shape, int)
        self.grid_A = np.asarray(grid_A, float)
        self.cell_A = np.asarray(cell_A, float)
        self.atoms = np.asarray(atoms, float)
        self.values = np.asarray(values, float)


def write_cube(cube, path):
    path = Path(path)
    scale = 1.0 / BOHR_TO_ANG
    with path.open("w") as f:
        f.write(cube.comment1 + "\n")
        f.write(cube.comment2 + "\n")
        o = cube.origin_A * scale
        f.write(f"{len(cube.atoms):5d} {o[0]:14.7f} {o[1]:14.7f} {o[2]:14.7f}\n")
        for n, v in zip(cube.shape, cube.grid_A):
            q = v * scale
            f.write(f"{int(n):5d} {q[0]:14.7f} {q[1]:14.7f} {q[2]:14.7f}\n")
        for row in cube.atoms:
            xyz = row[2:5] * scale
            f.write(f"{int(row[0]):5d} {row[1]:14.7f} {xyz[0]:14.7f} {xyz[1]:14.7f} {xyz[2]:14.7f}\n")
        flat = cube.values.ravel(order="C")
        for i in range(0, flat.size, 6):
            f.write(" ".join(f"{x:13.6e}" for x in flat[i:i+6]) + "\n")


def read_cube(path):
    path = Path(path)
    with path.open() as f:
        c1 = f.readline().rstrip()
        c2 = f.readline().rstrip()
        h = f.readline().split()
        nat = abs(int(h[0]))
        origin = np.array([float(x) for x in h[1:4]], float)
        shape = []
        vectors = []
        signs = []
        for _ in range(3):
            q = f.readline().split()
            n = int(q[0])
            signs.append(n)
            shape.append(abs(n))
            vectors.append([float(x) for x in q[1:4]])
        factor = BOHR_TO_ANG if all(x > 0 for x in signs) else 1.0
        origin *= factor
        grid = np.asarray(vectors, float) * factor
        atoms = []
        for _ in range(nat):
            q = f.readline().split()
            xyz = np.array([float(x) for x in q[2:5]], float) * factor
            atoms.append([int(float(q[0])), float(q[1]), *xyz])
        values = np.fromstring(f.read(), sep=" ")
    shape = np.asarray(shape, int)
    if values.size != int(np.prod(shape)):
        raise ValueError("Cube grid size does not match the number of values.")
    cell = grid * shape[:, None]
    return Cube(c1, c2, origin, shape, grid, cell, atoms, values.reshape(tuple(shape), order="C"))
