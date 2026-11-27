"""第五壳完整夸克见证可同时把六个强权行的两条重源投影严格归零。"""
from pathlib import Path
from runpy import run_path


NAME = "ch39-fifth-zero-strong-flavor"
FRAME = ((3,0,"i"),(3,1,"i"),(3,0,"r"),(2,0,"r"),(2,1,"r"),(2,0,"i"))
COEFFICIENTS = (299512305407,-42332840835,-45308455724,-291965817691,38888993493239,91781099086360,
                27519451906,-1603258322,7544390758,-23653065673,10607360171548,303634032291)
DIRECTIONS = ((14,1,"i"),(6,1,"r"),(7,1,"r"),(8,1,"r"),(7,1,"i"),(7,0,"r"),
              (8,1,"i"),(6,1,"i"),(14,1,"r"),(14,0,"r"),(6,4,"r"),(9,7,"r"),
              (8,4,"i"),(14,13,"i"),(6,0,"r"),(13,8,"i"),(13,7,"i"),(15,8,"r"),
              (14,10,"i"),(12,1,"r"),(9,6,"i"),(13,6,"i"),(17,12,"r"),(15,13,"r"),
              (9,1,"r"),(12,1,"i"),(9,8,"i"),(18,0,"r"),(18,9,"i"),(16,13,"i"),
              (14,2,"r"),(15,8,"i"),(18,7,"r"),(10,1,"r"),(12,0,"i"),(17,7,"r"))


def verify(root):
    def block_target(root, lib, blocks):
        return lib("fifth-free-frame-blocks")["solve"](root,lib,blocks,NAME,COEFFICIENTS,FRAME)

    chart=run_path(str(root/"ch39-fifth-flavor-cutoff-chart.py"))["solve"]
    return chart(root,NAME,(0,0,0,1),DIRECTIONS,frame=block_target,
                 rows=(0,9,15,16,17,18),target_radius="1e-13")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
