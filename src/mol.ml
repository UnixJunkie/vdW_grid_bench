(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module A = BatArray
module FSet = BatSet.Float
module Ht = BatHashtbl
module ISet = BatSet.Int
module L = BatList
module Log = Dolog.Log
module LO = Line_oriented

type t = { xs: float array; (* atoms X coords *)
           ys: float array; (* atoms Y coords *)
           zs: float array; (* atoms Z coords *)
           mutable center: V3.t;
           r_a: float array; (* vdW radii *)
           elt_a: int array; (* atomic numbers *)
           name: string } (* /!\ RENAMING MOLECULES IS _FORBIDDEN_ /!\ *)

let create mol_name num_atoms: t =
  let xs     = A.make num_atoms 0.0  in
  let ys     = A.make num_atoms 0.0  in
  let zs     = A.make num_atoms 0.0  in
  let center = V3.origin             in
  let r_a    = A.make num_atoms 0.0  in
  let elt_a  = A.make num_atoms 0    in
  { xs; ys; zs; center; r_a; elt_a; name = mol_name }

let num_atoms m: int =
  A.length m.r_a

let uget = A.unsafe_get

let set_radius m i r: unit =
  A.set m.r_a i r

let set_anum m i anum: unit =
  A.set m.elt_a i anum

let get_xs m: float array = m.xs
let get_ys m: float array = m.ys
let get_zs m: float array = m.zs

let get_xyz m i: V3.t =
  V3.make
    (uget m.xs i)
    (uget m.ys i)
    (uget m.zs i)

let get_all_atom_coords m: V3.t array =
  A.init (num_atoms m) (get_xyz m)

let get_anums m: int array =
  m.elt_a

let x_min_max m: float * float =
  A.min_max m.xs

let y_min_max m: float * float =
  A.min_max m.ys

let z_min_max m: float * float =
  A.min_max m.zs

(* mean_xyz *)
let get_center m: V3.t =
  m.center

let update_center m =
  m.center <- V3.{ x = A.favg m.xs ;
                   y = A.favg m.ys ;
                   z = A.favg m.zs }

(* only up to atom coordinates; a ligand's pqrs file will contain
   additional info *)
let pqrs_read_header_coords input =
  let header = input_line input in
  let num_atoms, num_rbonds, mol_name = Pqrs.parse_header header in
  (* Log.info "m:%s atoms:%d" mol_name num_atoms; *)
  let mol = create mol_name num_atoms in
  for i = 0 to num_atoms - 1 do
    let x, y, z, _q, r, elt = Pqrs.parse_body (input_line input) in
    set_radius mol i r;
    set_anum   mol i (Ptable.anum_of_symbol elt);
    A.set mol.xs i x;
    A.set mol.ys i y;
    A.set mol.zs i z
  done;
  update_center mol;
  (mol, num_rbonds)

let ligand_pqrs_read_one input =
  let mol, _num_rbonds = pqrs_read_header_coords input in
  mol

let first_ligand_of_pqrs_file (fn: string): t =
  LO.with_in_file fn ligand_pqrs_read_one
