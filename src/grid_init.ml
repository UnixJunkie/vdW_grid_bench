(* Copyright (C) 2024, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module A = Array
module BA = Bigarray
module BA3 = BA.Array3
module L = List

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

(* A generic vector: pull vector *)
type ('i,'a) vec = Vec of 'i *     (* exclusive upper bound for the index *)
                          ('i -> 'a)

(* An element-wise operation on a single vector *)
let vec_map : ('a -> 'b) -> ('i,'a) vec -> ('i,'b) vec =
 fun tf (Vec (n, f)) -> Vec (n, fun i -> tf (f i))
(* An element-wise binary vector operation *)
(* bounds should be physically equal
   (If they were not, one should have reconciled them first, or
   rolled into a element mapping)
*)
let zip_with : ('a -> 'b ->'c) -> ('i,'a) vec -> ('i,'b) vec -> ('i,'c) vec = 
  fun tf (Vec (n1,f1)) (Vec (n2,f2)) -> 
    assert( n1 == n2 );
    Vec (n1, fun i -> tf (f1 i) (f2 i))

type ('i,'a) reducer = 'a -> ('a -> 'a -> 'a) -> ('i,'a) vec -> 'a

let reducer_present : (int,'a) reducer = fun z f ->
  function Vec (upe,ix) ->
  let acc = ref z in 
    for i = 0 to upe-1 do
     acc := f !acc (ix i)
   done; !acc

let reducer_dyn : ('i code,'a code) reducer = fun z f -> 
  function Vec (upe,ix) ->
  .<let acc = ref .~z in 
    for i = 0 to .~upe-1 do
     acc := .~(f .<!acc>. (ix .<i>.))
   done; !acc>. 

(* Tabulating and memoizing the vector of parameters (protein number
   positions and LJ parameters with respect to the given ligand)
*)

let prot_vec_tabulate : (int,(float*float*float)*(float*float)) vec ->
  (int * Float.Array.t) = function Vec (n,ix) ->
    List.init n (fun i ->
      let ((x,y,z),(xi,di)) = ix i in
      Float.Array.of_list [x;y;z;xi;di]) 
    |> Float.Array.concat 
    |> fun a -> (n,a)

let to_prot_vec : (int * Float.Array.t) -> 
  (int,(float*float*float)*(float*float)) vec = fun (n,arr) ->
    Vec (n, fun i' ->
      let i = 5*i' in
      let x = Float.Array.unsafe_get arr i in
      let y = Float.Array.unsafe_get arr (i+1) in
      let z = Float.Array.unsafe_get arr (i+2) in
      let xi = Float.Array.unsafe_get arr (i+3) in
      let di = Float.Array.unsafe_get arr (i+4) in
      ((x,y,z),(xi,di)))

let prot_vec_memo : (int,(float*float*float)*(float*float)) vec ->
   (int,(float*float*float)*(float*float)) vec =
  fun v -> v |> prot_vec_tabulate |> to_prot_vec

let to_prot_vec_staged : (int code * Float.Array.t code) -> 
  (int code,(float code*float code*float code)*(float code*float code)) vec =
  fun (n,arr) ->
    let arr_get i = .<Float.Array.unsafe_get .~arr .~i>. in
    Vec (n, fun i' ->
      let i = genlet .<5* .~i'>. in
      let x = arr_get i in
      let y = arr_get  .<.~i+1>. in
      let z = arr_get  .<.~i+2>. in
      let xi = arr_get .<.~i+3>. in
      let di = arr_get .<.~i+4>. in
      ((x,y,z),(xi,di)))

(* Make staged code easier to read  *)
module C = struct
  let zero = .<0.>. 
  let ( +. ) = fun x y -> .<.~x +. .~y>.
  let ( -. ) = fun x y -> .<.~x -. .~y>.
  let ( *. ) = fun x y -> .<.~x *. .~y>.
  let ( /. ) = fun x y -> .<.~x /. .~y>.
  let sqr x = letl x @@ fun x -> x *. x
  let sqrt x = .<sqrt .~x>.
  let pow6 x =
   letl x @@ fun x -> letl (x *. x) @@ fun x -> 
   genlet (x *. x *. x)
  let vec_get : 'a array code -> int code -> 'a code = fun v i ->
   .< Array.unsafe_get .~v .~i >.              (* or unsafe_get *)
  let vec_set : 'a array code -> int code -> 'a code -> unit code = fun v i x ->
   .< Array.unsafe_set .~v .~i .~x >.          (* or unsafe_set *)
end

let v3_make_staged : float code -> float code -> float code ->
   (float code * float code * float code) =
   fun x y z -> (x,y,z)

let v3_dist : (float code * float code * float code) -> 
              (float code * float code * float code) -> float code = 
   fun (v11,v12,v13) (v21,v22,v23) ->
    C.(sqrt (sqr (v11 -. v21) +. sqr (v12 -. v22) +. sqr (v13 -. v23)))
let non_zero_dist_staged : float code -> float code = fun x ->
    letl x @@ fun x -> .<if .~x < 0.01 then 0.01 else .~x >.


let timeit : string -> (unit -> 'a) -> 'a = fun msg thunk ->
   let tbeg = Sys.time () in
   let r = thunk () in
   Printf.printf "%s took %g sec\n" msg (Sys.time () -. tbeg); 
   r

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


  (* To confirm that different versions of code in fact produce the same
     grid (or, at least, seemingly the same to some degree of confidence)
  *)
  let fingerprint {grid;nx;ny;nz} : float =
    let size = nx * ny * nz in
    let m1 = BA.(reshape_1 (genarray_of_array3 grid) size) in
    let sum = ref 0.0 in
    for i = 0 to size -1 do sum := !sum +. BA.Array1.get m1 i done;
    !sum

  (* The real code
  let to_file (fn: string) (x: t): unit =
    let c = open_out fn in
    match Marshal.to_channel c x [Marshal.No_sharing] with
    | exception e -> close_out c; raise e
    | _ -> close_out c
 *)
  (* Used for benchmarking *)
  let to_file (fn: string) (x: t): unit =
    Printf.printf "Prepared grid: fingerprint %.17g\n" (fingerprint x)

   (*
  let from_file (fn: string): t =
    LO.restore fn
  *)

  (* initialize g.grid as a vdW interaction grid for ligand atomic number [l_a]
     WARNING: this code needs to go fast, because the grid might be big
              if discretization step is small and/or protein is big *)
  let vdW_grid (nx,ny,nz) (x_min, y_min, z_min) dx prot_vec =
    let prot_vec = prot_vec_memo prot_vec in
    let grid = BA3.create BA.float32 BA.c_layout nx ny nz in
    let size = nx * ny * nz in
    let m1 = BA.(reshape_1 (genarray_of_array3 grid) size) in
    let ix = ref 0 in
    let x = ref x_min in
    for i = 0 to nx - 1 do
      let y = ref y_min in
      for j = 0 to ny - 1 do
        let z = ref z_min in
        for k = 0 to nz - 1 do
          let l_p = V3.make !x !y !z in (* ligand atom position *)
          vec_map (fun ((x,y,z),(x_ij,d_ij)) ->
            let r_ij = non_zero_dist (V3.dist l_p (V3.make x y z)) in
            let p6 = pow6 (x_ij /. r_ij) in
            (d_ij *. ((-2.0 *. p6) +. (p6 *. p6)))) prot_vec
           |> reducer_present 0.0 (+.) |> BA.Array1.unsafe_set m1 !ix;
        incr ix;
        z := !z +. dx
        done; y := !y +. dx
      done; x := !x +. dx
    done;
    grid

  (* specialized to the grid parameters. The prot vector is dynamic *)
  (* it is deliberately made to be very similar to vdW_grid *)
  let vdW_grid_staged ((nx,ny,nz):(int*int*int)) 
   ((x_min, y_min, z_min):(float*float*float)) (dx:float) prot_vec =
   .<
    let grid = Bigarray.(Array3.create float32 c_layout nx ny nz) in
    let size = nx * ny * nz in
    let m1 = Bigarray.(reshape_1 (genarray_of_array3 grid) size) in
    let ix = ref 0 in
    let x = ref x_min in
    for i = 0 to nx - 1 do
      let y = ref y_min in
      for j = 0 to ny - 1 do
        let z = ref z_min in
        for k = 0 to nz - 1 do
        .~(
          let l_p = v3_make_staged .<!x>. .<!y>. .<!z>. in (* ligand at pos *)
          vec_map (fun (p_p,(x_ij,d_ij)) ->
            let open C in
            let r_ij = non_zero_dist_staged (v3_dist l_p p_p) in
            let p6 = genlet @@ pow6 (x_ij /. r_ij) in
            d_ij *. ((.<-2.0>. *. p6) +. (p6 *. p6))) prot_vec
           |> reducer_dyn C.zero C.(+.)
           |> fun v -> .<Bigarray.Array1.unsafe_set m1 !ix .~v>.);
        incr ix;
        z := !z +. dx
        done; y := !y +. dx
      done; x := !x +. dx
    done;
    grid
   >.

  (* run-time code generation *)
  let vdW_grid_gen :
      (int*int*int) -> (float*float*float) -> float -> (int,_) vec ->
      _ BA3.t =
    let memo = ref [] in
    fun ns mins dx prot_vec ->
      let prot = prot_vec_tabulate prot_vec in
      match List.assoc_opt (ns,mins) !memo with
      | Some f -> 
          Printf.printf "Reusing the memoized grid filling code\n";
          f prot
      | None ->
          let cde = 
            timeit "generating code" @@ fun () ->
            .<fun (n,prot) ->
              .~(vdW_grid_staged ns mins dx 
                   (to_prot_vec_staged (.<n>.,.<prot>.)))>.
            in
          print_code Format.std_formatter cde;
          let f = timeit "compiling generated code" @@ fun () ->
            Runnative.run cde
          in
          memo := ((ns,mins),f) :: !memo;
          f prot

  let create step low high prot_vec =
    let x_min, y_min, z_min = V3.to_triplet low  in
    let x_max, y_max, z_max = V3.to_triplet high in
    let dx = x_max -. x_min in
    let dy = y_max -. y_min in
    let dz = z_max -. z_min in
    let nx = int_of_float (ceil (dx /. step)) in
    let ny = int_of_float (ceil (dy /. step)) in
    let nz = int_of_float (ceil (dz /. step)) in
    Printf.printf "Grid.create: nx,ny,nz=%d,%d,%d\n" nx ny nz;
    { low; high; step; nx; ny; nz; 
      grid=
      vdW_grid (nx,ny,nz) (x_min,y_min,z_min) step prot_vec}

  let create_gen step low high prot_vec =
    let x_min, y_min, z_min = V3.to_triplet low  in
    let x_max, y_max, z_max = V3.to_triplet high in
    let dx = x_max -. x_min in
    let dy = y_max -. y_min in
    let dz = z_max -. z_min in
    let nx = int_of_float (ceil (dx /. step)) in
    let ny = int_of_float (ceil (dy /. step)) in
    let nz = int_of_float (ceil (dz /. step)) in
    Printf.printf "Grid.create: nx,ny,nz=%d,%d,%d\n" nx ny nz;
    { low; high; step; nx; ny; nz; 
      grid=
      vdW_grid_gen (nx,ny,nz) (x_min,y_min,z_min) step prot_vec}
end

(* parse pqrs file only made of atom lines *)
let read_prot_pqrs_file_no_header fn =
  let lines = 
    let c = open_in fn in
    let rec loop acc = 
      match In_channel.input_line c with
      | Some s -> loop (s::acc) 
      | None -> List.rev acc
    in match loop [] with
    | exception e -> close_in c; raise e
    | l -> close_in c; l
  in
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
  let _verbose = ref false in
  let lig_fn = ref "" in
  let prot_fn = ref "" in
  let step = ref 0.0 in
  let _nprocs = ref 1 in
  let speclist = 
    [("-verbose", Arg.Set _verbose, "Output debug information");
      ("-l", Arg.Set_string lig_fn, "ligand <FILE.pqrs>");
      ("-p", Arg.Set_string prot_fn, "protein <FILE.pqrs>");
      ("-dx", Arg.Set_float step, "grid step");
      ("-np", Arg.Set_int _nprocs, "nprocs (default=1)");] in
  let () = Arg.parse speclist (fun _ -> failwith "unexpected extra args")
    "usage:\n  \
              %s\n  \
              -l <FILE.pqrs>: ligand\n  \
              -p <FILE.pqrs>: protein\n  \
              -dx <float>: grid step\n  \
              [-np <int>]: nprocs (default=1)\n  \
              [-v]: verbose/debug mode\n"
  in
  let lig_fn = !lig_fn in
  let prot_fn = !prot_fn in
  let step = 
    if !step = 0.0 then 0.2 else
    (assert(!step > 0.01); (* because of non_zero_dist *)
     !step)
  in
  let _nprocs = !_nprocs in
  (* END CLI parsing ------------------------------------------------------- *)
  let lig = Mol.first_ligand_of_pqrs_file lig_fn in
  let prot = read_prot_pqrs_file_no_header prot_fn in
  let x_min, x_max = Mol.x_min_max prot in
  let y_min, y_max = Mol.y_min_max prot in
  let z_min, z_max = Mol.z_min_max prot in
  let low_corner = V3.make x_min y_min z_min in
  let high_corner = V3.make x_max y_max z_max in
  Printf.printf "%d atoms in %s\n" (Mol.num_atoms lig)  lig_fn;
  Printf.printf "%d atoms in %s\n" (Mol.num_atoms prot) prot_fn;
  let lig_anums =
    let anums = A.to_list (Mol.get_anums lig) in
    L.sort_uniq Int.compare anums in
  Printf.printf "%d anums in %s\n" (L.length lig_anums) lig_fn;
  let prot_coords = Mol.get_all_atom_coords prot in
  let prot_anums = Mol.get_anums prot in  
  (* initialize all vdW grids *)
  let vdW_grids =
    L.map (fun l_a ->
        Printf.printf "l_a: %d\n" l_a;
        let prot_vec = 
          Vec (Mol.num_atoms prot,
               fun i -> ((prot_coords.(i) |> V3.to_triplet),
                         let vdw = UFF.vdW_xiDi l_a prot_anums.(i) in
                         (vdw.x_ij,vdw.d_ij)))
        in
        (l_a,
        Grid.create_gen step low_corner high_corner prot_vec)
      ) lig_anums in
  (* dump to disk *)
  L.iter (fun (l_a, g) ->
      let fn = sprintf "%s.%d\n" lig_fn l_a in
      Printf.printf "creating %s" fn;
      Grid.to_file fn g
    ) vdW_grids

let () = main ()
