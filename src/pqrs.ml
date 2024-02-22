(* Copyright (C) 2023, Francois Berenger
 * Tsuda laboratory, The University of Tokyo,
 * 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan. *)

module Log = Dolog.Log

let parse_header line: int * int * string =
  try (* ligand header w/ rot-bonds *)
    Scanf.sscanf line "%d:%d:%s" (fun num_atoms num_bonds name ->
        (num_atoms, num_bonds, name)
      )
  with _exn ->
    begin
      try (* receptor header (rot-bonds are not mentioned) *)
        Scanf.sscanf line "%d:%s" (fun n name -> (n, -1, name))
      with exn ->
        let () = Log.error "Mol.parse_pqrs_header: cannot parse: %s" line in
        raise exn
    end

let parse_body line: float * float * float * float * float * string =
  try Scanf.sscanf line "%f %f %f %f %f %s"
        (fun x y z q r elt -> (x, y, z, q, r, elt))
  with exn ->
    let () =
      Log.error "Molecule.parse_pqrs_body: cannot parse: %s" line in
    raise exn
