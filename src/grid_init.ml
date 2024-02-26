(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module A = BatArray
module BA = BatBigarray
module BA3 = BA.Array3
module CLI = Minicli.CLI
module L = BatList
module LO = Line_oriented
module Log = Dolog.Log

open Printf

let pow3 x =
  x *. x *. x

let pow6 x =
  pow3 (x *. x)

(* avoid division by 0 in case of too small inter atomic distance *)
let non_zero_dist x =
  if x < 0.01 then
    0.01 (* trick used in D. Horvath's S4MPLE docking program *)
  else
    x

let non_zero_dist2 x =
  if x < 0.0001 then
    0.0001
  else
    x

module Grid = struct

  type t = {
    (* meta data *)
    low: V3.t;
    high: V3.t;
    step: float; (* Cubic grid: dx = dy = dz *)
    (* integer grid dimensions *)
    nx: int;
    ny: int;
    nz: int;
    (* data *)
    grid: (float, BA.float32_elt, BA.c_layout) BA3.t }

  let create step low high =
    let x_min, y_min, z_min = V3.to_triplet low  in
    let x_max, y_max, z_max = V3.to_triplet high in
    let dx = x_max -. x_min in
    let dy = y_max -. y_min in
    let dz = z_max -. z_min in
    let nx = int_of_float (ceil (dx /. step)) in
    let ny = int_of_float (ceil (dy /. step)) in
    let nz = int_of_float (ceil (dz /. step)) in
    Log.info "Grid.create: nx,ny,nz=%d,%d,%d" nx ny nz;
    let grid = BA3.create BA.float32 BA.c_layout nx ny nz in
    (* init w/ 0s *)
    BA3.fill grid 0.0;
    { low; high; step; nx; ny; nz; grid }

  let to_file (fn: string) (x: t): unit =
    LO.save fn x

  let from_file (fn: string): t =
    LO.restore fn

  (* initialize g.grid as a vdW interaction grid for ligand atomic number [l_a]
     WARNING: this code needs to go fast, because the grid might be big
              if discretization step is small and/or protein is big *)
  let vdW_grid g prot_atoms l_a =
    let x_min, y_min, z_min = V3.to_triplet g.low in
    let dx = g.step in
    let n = A.length prot_atoms in
    for i = 0 to g.nx - 1 do
      let x = x_min +. (float i) *. dx in
      for j = 0 to g.ny - 1 do
        let y = y_min +. (float j) *. dx in
        for k = 0 to g.nz - 1 do
          let z = z_min +. (float k) *. dx in
          (* ligand atom position *)
          let l_p = V3.make x y z in
          let grid_i_j_k = ref (BA3.unsafe_get g.grid i j k) in
          for m = 0 to n - 1 do (* over all protein atoms *)
            let p_p, p_a = A.unsafe_get prot_atoms m in
            let r_ij2 = non_zero_dist2 (V3.dist2 l_p p_p) in
            let vdw = UFF.vdW_xiDi l_a p_a in
            let p6 = pow3 (vdw.x_ij *. vdw.x_ij /. r_ij2) in
            (* g.grid.{i,j,k} <-
               g.grid.{i,j,k} +.
               (vdw.d_ij *. ((-2.0 *. p6) +. (p6 *. p6))) *)
            grid_i_j_k :=
              (!grid_i_j_k +.
               (vdw.d_ij *. (p6 *. (p6 -. 2.0))))
          done;
          BA3.unsafe_set g.grid i j k !grid_i_j_k
        done
      done
    done;
    g (* for more functional programming style *)

end

(* parse pqrs file only made of atom lines *)
let read_prot_pqrs_file_no_header fn =
  let lines = LO.lines_of_file fn in
  let num_atoms = L.length lines in
  let mol = Mol.create fn num_atoms in
  L.iteri (fun i line ->
      let x, y, z, _q, r, elt = Pqrs.parse_body line in
      Mol.set_radius mol i r;
      Mol.set_anum   mol i (Ptable.anum_of_symbol elt);
      A.set mol.xs i x;
      A.set mol.ys i y;
      A.set mol.zs i z
    ) lines;
  Mol.update_center mol;
  mol

let main () =
  Log.(set_log_level INFO);
  Log.color_on ();
  Log.(set_prefix_builder short_prefix_builder);
  let argc, args = CLI.init () in
  if argc = 1 then
    (eprintf "usage:\n  \
              %s\n  \
              -l <FILE.pqrs>: ligand\n  \
              -p <FILE.pqrs>: protein\n  \
              -dx <float>: grid step\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-v]: verbose/debug mode\n"
       Sys.argv.(0);
     exit 1);
  let _verbose = CLI.get_set_bool ["-v"] args in
  let lig_fn = CLI.get_string ["-l"] args in
  let prot_fn = CLI.get_string ["-p"] args in
  let step = CLI.get_float_def ["-dx"] args 0.2 in
  assert(step > 0.01); (* because of non_zero_dist *)
  let _nprocs = CLI.get_int_def ["-np"] args 1 in
  CLI.finalize ();
  (* END CLI parsing ------------------------------------------------------- *)
  let lig = Mol.first_ligand_of_pqrs_file lig_fn in
  let prot = read_prot_pqrs_file_no_header prot_fn in
  let x_min, x_max = Mol.x_min_max prot in
  let y_min, y_max = Mol.y_min_max prot in
  let z_min, z_max = Mol.z_min_max prot in
  let low_corner = V3.make x_min y_min z_min in
  let high_corner = V3.make x_max y_max z_max in
  Log.info "%d atoms in %s" (Mol.num_atoms lig)  lig_fn;
  Log.info "%d atoms in %s" (Mol.num_atoms prot) prot_fn;
  let lig_anums =
    let anums = A.to_list (Mol.get_anums lig) in
    L.sort_uniq BatInt.compare anums in
  Log.info "%d anums in %s" (L.length lig_anums) lig_fn;
  let prot_coords = Mol.get_all_atom_coords prot in
  let prot_anums = Mol.get_anums prot in
  let vdW_prot_atoms = A.combine prot_coords prot_anums in
  (* initialize all vdW grids *)
  let vdW_grids =
    L.map (fun l_a ->
        Log.info "l_a: %d" l_a;
        let g = Grid.create step low_corner high_corner in
        (l_a, Grid.vdW_grid g vdW_prot_atoms l_a)
      ) lig_anums in
  (* dump to disk *)
  L.iter (fun (l_a, g) ->
      let fn = sprintf "%s.%d" lig_fn l_a in
      Log.info "creating %s" fn;
      Grid.to_file fn g
    ) vdW_grids

let () = main ()
