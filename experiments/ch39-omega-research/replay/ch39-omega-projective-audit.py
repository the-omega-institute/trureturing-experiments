"""Rebuild the pinned Omega projective-circle imports and print their axiom sets."""
import argparse
import hashlib
import os
import re
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

REVISION = "778fa34bb890f75a63d5e2791fe3e1c56707a9d6"
THEOREMS = ("kappa_scale", "kappa_add_spread_eq_one", "cayley_circle_unit", "cayley_compose")


def audit(checkout):
    src = checkout / "src"
    built = set()
    with TemporaryDirectory() as directory:
        build = Path(directory).resolve()
        env = {**os.environ, "LEAN_PATH": str(build)}

        def compile_module(module):
            if module in built:
                return
            path = Path(*module.split(".")).with_suffix(".lean")
            source = (src / path).read_bytes()
            digest = hashlib.sha1(b"blob " + str(len(source)).encode() + b"\0" + source).hexdigest()
            expected = subprocess.check_output(
                ["git", "rev-parse", f"{REVISION}:src/{path.as_posix()}"], cwd=checkout, text=True).strip()
            assert expected == digest
            for dep in re.findall(r"^(?:public |private )?import (Omega\w*(?:\.\w+)*)", source.decode(), re.M):
                compile_module(dep)
            obj = (build / path).with_suffix(".olean")
            obj.parent.mkdir(exist_ok=True, parents=True)
            subprocess.run(["lean", "-R", str(src), "-o", str(obj), str(src/path)],
                           cwd=checkout, env=env, check=True)
            built.add(module)

        compile_module("OmegaMath.Projective.Circle")
        statements = "".join("#print axioms Omega.Projective.Circle." + x + "\n" for x in THEOREMS)
        result = subprocess.check_output(["lean", "--stdin"], cwd=checkout, env=env, text=True,
                                         input="prelude\nimport OmegaMath.Projective.Circle\n" + statements)
        assert all(f"'Omega.Projective.Circle.{name}' does not depend on any axioms" in result for name in THEOREMS)
        print(result, end="")
        print(len(built), "pinned source modules rebuilt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("checkout", type=Path, help="trureturing checkout with its pinned Lean toolchain")
    audit(parser.parse_args().checkout.resolve())
