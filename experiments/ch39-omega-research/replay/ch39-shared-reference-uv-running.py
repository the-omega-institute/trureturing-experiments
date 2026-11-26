"""Run the same gauge thresholds with only the minimum vertex-response states."""
from pathlib import Path
from runpy import run_path


if __name__=="__main__":
    run=run_path(str(Path(__file__).with_name("ch39-kinetic-rh-uv-running.py")))
    run["main"](shared=True)
