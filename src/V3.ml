(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

(* inserted the needed parts from Vector3 *)
type t = {x: float; y:float; z:float}

let make x y z =
  { x = x ;
    y = y ;
    z = z }

let origin =
  { x = 0. ;
    y = 0. ;
    z = 0. }

let to_triplet v =
  (v.x, v.y, v.z)

(* for a .xyz trajectory file *)
let xyz_dump out v =
  Printf.fprintf out " %.3f %.3f %.3f\n" v.x v.y v.z

(* [dist2]: squared distance,
   saves one call to sqrt compared to [dist] below *)
let dist2 u v =
  (*    dx2 *)
  (u.x -. v.x) *.
  (u.x -. v.x)
  +. (* dy2 *)
  (u.y -. v.y) *.
  (u.y -. v.y)
  +. (* dz2 *)
  (u.z -. v.z) *.
  (u.z -. v.z)

(* the one in library vector3 should be optimized *)
let dist u v =
  sqrt (dist2 u v)
