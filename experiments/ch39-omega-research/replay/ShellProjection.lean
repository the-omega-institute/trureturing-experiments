import Init
set_option maxRecDepth 80000
set_option maxHeartbeats 4000000

structure H where
  w : Int
  x : Int
  y : Int
  z : Int

def eqInt (a b : Int) : Bool :=
  Int.casesOn a
    (fun n => Int.casesOn b (fun m => Nat.beq n m) (fun _ => false))
    (fun n => Int.casesOn b (fun _ => false) (fun m => Nat.beq n m))

def same (a b : H) : Bool :=
  eqInt a.w b.w && eqInt a.x b.x && eqInt a.y b.y && eqInt a.z b.z

def checkAll (f : H → Bool) : List H → Bool
  | [] => true
  | a :: rest => f a && checkAll f rest

def plus (a b : H) : H := ⟨a.w+b.w,a.x+b.x,a.y+b.y,a.z+b.z⟩
def smul (n : Int) (a : H) : H := ⟨n*a.w,n*a.x,n*a.y,n*a.z⟩
def mul (a b : H) : H :=
  ⟨a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z, a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,
   a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x, a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w⟩
def dot (a b : H) : Int := a.x*b.x+a.y*b.y+a.z*b.z
def axis : H := ⟨0,0,1,4⟩
def plane (p v : H) : H := plus (smul 17 v) (mul (mul p v) axis)
def basis : List H := [⟨1,0,0,0⟩,⟨0,1,0,0⟩,⟨0,0,1,0⟩,⟨0,0,0,1⟩]
def shell : List H := [⟨0,0,1,4⟩,⟨0,0,4,1⟩,⟨0,1,0,4⟩,⟨0,1,4,0⟩,⟨0,4,0,1⟩,⟨0,4,1,0⟩,⟨0,2,2,3⟩,⟨0,2,3,2⟩,⟨0,3,2,2⟩]
def compress (p q : H) : Bool := checkAll (fun v =>
  same (plane p (plane q (plane p v))) (smul (34*(17+dot p q)) (plane p v))) basis
def oriented : Bool :=
  let p : H := ⟨0,3,2,2⟩
  let q : H := ⟨0,2,3,2⟩
  let r : H := ⟨0,2,2,3⟩
  checkAll (fun v =>
    let z := plane p v
    same (plane p (plane q (plane r z))) (smul 34 (plus (smul 1105 z) (smul 7 (mul z axis))))) basis
theorem shell_compressions : (checkAll (fun p => checkAll (compress p) shell) shell) = true := rfl
theorem inner_oriented_triangle : oriented = true := rfl
#print axioms shell_compressions
#print axioms inner_oriented_triangle

def innerAxis : Bool :=
  let p : H := ⟨0,3,2,2⟩
  let q : H := ⟨0,2,3,2⟩
  let r : H := ⟨0,2,2,3⟩
  let u := fun (a x : H) => mul (mul a x) axis
  let triple := fun (a b c x : H) => u a (u b (u c x))
  checkAll (fun x =>
    let forward := plus (triple p q r x) (plus (triple q r p x) (triple r p q x))
    let reverse := plus (triple p r q x) (plus (triple r q p x) (triple q p r x))
    same (plus forward (smul (-1) reverse)) (smul 714 (mul x axis))) basis

theorem oriented_axis_from_shell : innerAxis = true := rfl
#print axioms oriented_axis_from_shell
