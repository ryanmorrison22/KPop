(*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*)

open BiOCamLib
open Better
open KPop

module RegisterType =
  struct
    type t =
      | Twister
      | Twisted
    let of_string = function
      | "T" -> Twister
      | "t" -> Twisted
      | s ->
        Printf.sprintf "(%s): Invalid register type '%s'" __FUNCTION__ s |> failwith
  end

module KeepAtMost =
  struct
    type t = int option
    let of_string = function
      | "all" -> None
      | s ->
        try
          let res = int_of_string s in
          if res <= 0 then
            raise Exit;
          Some res
        with _ ->
          Printf.sprintf "(%s): Invalid keep_at_most '%s'" __FUNCTION__ s |> failwith
    let to_string = function
      | None -> "all"
      | Some n -> string_of_int n
  end

type to_do_t =
  | Empty of RegisterType.t
  | Binary_to_register of RegisterType.t * string (* Input prefix *)
  | Tables_to_register of RegisterType.t * string (* Input prefix *)
  | Add_binary_to_twisted of string (* Input prefix *)
  | Twist_database of string (* Input prefix *)
(* | Add_kmers_binary_to_twisted of string *)
  | Register_to_binary of RegisterType.t * string (* Output prefix *)
  | Set_precision_tables of int
  | Set_precision_splits of int
  | Register_to_tables of RegisterType.t * string (* Output prefix *)
  | Set_metric of Space.Distance.Metric.t
  | Set_distance of Space.Distance.t
  | Set_distance_normalize of bool
  (* Computes embeddings from twisted vectors using the current metric/distance.
     Parameter is output prefix *)
  | Embeddings_from_twisted of string
  | Set_splits_algorithm of Twisted.SplitsAlgorithm.t
  | Set_splits_keep_at_most of int
  | Splits_from_twisted of string (* Output prefix *)
  | Set_summary_keep_at_most of KeepAtMost.t
  (* Parameters are input probes, output prefix for summary/distance matrix,
      and a boolean specifying whether the distance matrix should be output *)
  | Summary_from_twisted_binary of string * string * bool
  | Set_neighbors_keep_at_most of KeepAtMost.t
  | Set_neighbors_guard_policy of Twisted.NeighborsPolicy.t
  | Set_neighbors_index_type of Interfaiss.Type.t
  (* Parameters are input probes and output summary *)
  | Summary_from_twisted_neighbors of string * string

module Defaults =
  struct
    let distance = Space.Distance.of_string "euclidean"
    let distance_normalize = false
    let metric = Space.Distance.Metric.of_string "powers(1,1,1)"
    let precision_tables = 15
    let precision_splits = 10
    let splits_algorithm = Twisted.SplitsAlgorithm.of_string "gaps"
    let splits_keep_at_most = 10000
    let summary_keep_at_most = Some 2
    let neighbors_keep_at_most = Some 6
    let neighbors_guard_policy = Twisted.NeighborsPolicy.of_string "times(2)"
    let neighbors_index_type = Interfaiss.Type.of_string "hnsw(32)"
  end

module Parameters =
  struct
    let program = ref []
    let threads = Processes.Parallel.get_nproc () |> ref
    let debug_twisting = ref false
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopTwistDB";
  version = "46";
  date = "03-Nov-2025"
} and authors = [
  "2022-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com";
  "2024     ", "Ünsal Öztürk", "uensal.oeztuerk@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    TA.make_separator_multiline [ ""; "Actions on the database registers - Input/Output operations:" ];
    [ "-0"; "--zero"; "--empty" ],
      Some "'T'|'t'",
      [ "load an empty database into the specified register";
        " ('T'=twister; 't'=twisted)" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Twisted as register_type ->
          Empty register_type |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "'T'|'t' <binary_file_prefix>",
      [ "load the specified binary database into the specified register";
        " ('T'=twister; 't'=twisted).";
        "File extension is automatically determined depending on database type";
        " (will be '.KPopTwister' or '.KPopTwisted', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Twisted as register_type ->
          Binary_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-I"; "--Input" ],
      Some "'T'|'t' <tabular_file_prefix>",
      [ "load the specified tabular database(s) into the specified register";
        " ('T'=twister; 't'=twisted).";
        "File extension is automatically determined depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt;';
                or: '.KPopInertia.txt' and '.KPopTwisted.txt', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Twisted as register_type ->
          Tables_to_register (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-a"; "--add"; "--add-to-twisted" ],
      Some "<binary_file_prefix>",
      [ "add the contents of the specified binary database to the twisted register.";
        "File extension is automatically determined";
        " (will be '.KPopTwisted', unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Add_binary_to_twisted (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "'T'|'t' <binary_file_prefix>",
      [ "save the database present in the specified register";
        " ('T'=twister; 't'=twisted)";
        "to the specified binary file.";
        "File extension is automatically assigned depending on database type";
        " (will be '.KPopTwister' or '.KPopTwisted', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        match TA.get_parameter () |> RegisterType.of_string with
        | Twister | Twisted as register_type ->
          Register_to_binary (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--precision-for-tables" ],
      Some "<positive_integer>",
      [ "set how many precision digits should be used when outputting numbers";
        "in tabular formats" ],
      TA.Default (fun () -> string_of_int Defaults.precision_tables),
      (fun _ -> Set_precision_tables (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-O"; "--Output" ],
      Some "'T'|'t' <tabular_file_prefix>",
      [ "save the database present in the specified register";
        " ('T'=twister; 't'=twisted)";
        "to the specified tabular files.";
        "File extensions are automatically assigned depending on database type";
        " (will be: '.KPopTwister.txt' and '.KPopInertia.txt';";
        "       or: '.KPopInertia.txt' and '.KPopTwisted.txt', respectively,";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let register_type = TA.get_parameter () |> RegisterType.of_string in
        Register_to_tables (register_type, TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Actions on the database register - Other operations:" ];
    [ "-t"; "--twist"; "--twist-kmers"; "--twist-spectra" ],
      Some "<binary_file_prefix>",
      [ "twist the k-mer spectra contained in the specified binary database";
        "according to the transformation present in the twister register,";
        "and add the results to the database loaded in the twisted register" ],
      TA.Optional,
      (fun _ -> Twist_database (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-m"; "--metric"; "--metric-function" ],
      Some "'flat'|'powers('POWERS_PARAMETERS')'",
      [ "where POWERS_PARAMETERS := ";
        "        INTERNAL_POWER','FRACTIONAL_ACCUMULATIVE_THRESHOLD','EXTERNAL_POWER";
        "      INTERNAL_POWER := <non-negative_float>";
        "      FRACTIONAL_ACCUMULATIVE_THRESHOLD := <fractional_float>";
        "      EXTERNAL_POWER := <non-negative_float>";
        "Set the metric function to be used when computing distances.";
        "The 'power' transformation is computed as follows:";
        " (1) the inertia vector is raised to INTERNAL_POWER and normalized;";
        " (2) elements are summed in order until FRACTIONAL_ACCUMULATIVE_THRESHOLD";
        "     (a number between 0. and 1.) is reached, while the elements";
        "     above the threshold are set to zero";
        " (3) the resulting vector is raised to EXTERNAL_POWER and normalized.";
        "Note that";
        " 'flat'";
        "(which is equivalent to 'power(0,1,1)' or 'power(1,1,0)')";
        "disregards inertia, i.e. it is the same as standard coordinates, while";
        " 'power(1,1,1)'";
        "leaves inertia unchanged, i.e. it is the same as principal coordinates" ],
      TA.Default (fun () -> Space.Distance.Metric.to_string Defaults.metric),
      (fun _ ->
        Set_metric (TA.get_parameter () |> Space.Distance.Metric.of_string) |> List.accum Parameters.program);
    [ "--distance"; "--distance-function" ],
      Some "'euclidean'|'cosine'|'angle'|'minkowski('<non-negative_float>')'",
      [ "set the function to be used when computing distances.";
        "The parameter for 'minkowski' is the power.";
        "Note that:";
        " 'euclidean' is the same as 'minkowski(2)';";
        " 'cosine' is the same as ('euclidean'^2)/2, or 1 - cos theta;";
        " 'angle' is the same as arccos(1 - ('euclidean'^2)/2), or theta,";
        "where theta is the relative angle between the two embeddings" ],
      TA.Default (fun () -> Space.Distance.to_string Defaults.distance),
      (fun _ -> Set_distance (TA.get_parameter () |> Space.Distance.of_string) |> List.accum Parameters.program);
    [ "--distance-normalize"; "--distance-normalization" ],
      Some "'true'|'false'",
      [ "whether to normalize twisted vectors before computing distances.";
        "It must be 'true' when the distance function is 'cosine' or 'angle'" ],
      TA.Default (fun () -> string_of_bool Defaults.distance_normalize),
      (fun _ -> Set_distance_normalize (TA.get_parameter_boolean ()) |> List.accum Parameters.program);
    [ "-e"; "--embeddings"; "--compute-embeddings"; "--twisted-to-embeddings" ],
      Some "<tabular_file_prefix>",
      [ "compute embeddings from the vectors present in the twisted register";
        "using the current metric function, distance function and normalization.";
        "The result will be written to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopVectors' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Embeddings_from_twisted (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "--distances-summarize-at-most"; "--distances-in-summary" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest sequences to be printed";
        "when summarizing distances.";
        "Note that more might be printed anyway in case of ties.";
        "The statistics in the summary will be computed on all sequences" ],
      TA.Default (fun () -> KeepAtMost.to_string Defaults.summary_keep_at_most),
      (fun _ ->
        Set_summary_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
    [ "-d"; "--summarize-distances"; "--compute-and-summarize-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "for each vector present in the twisted register, compute distances";
        "to all vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the current metric function, distance function and normalization;";
        "summarize them and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter (), false) |> List.accum Parameters.program);
    [ "-D"; "--summarize-and-output-distances"; "--compute-summarize-and-output-distances" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "same as option '-d', but additionally output the distance matrix";
        "in tabular form.";
        "File extensions are automatically assigned";
        " (will be '.KPopSummary.txt' and '.KPopDMatrix.txt',";
        "  unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_binary (twisted_prefix, TA.get_parameter (), true) |> List.accum Parameters.program);
    [ "--neighbors-index-type"; "--neighbors-faiss-index-type" ],
      Some "'flat'|'pq('PQ_PARAMETERS')'|'hnsw('<positive_integer>')'",
      [ "where PQ_PARAMETERS :=";
        " <positive_integer>','<positive_integer>' :";
        "set the type of Faiss index used to compute nearest neighbors.";
        "Parameters for 'pq' are:";
        " number of subquantizers; bits per subquantizer.";
        "Note that the product of the two must be less than or equal to";
        "the number of dimensions of the twisted vectors.";
        "The parameter for 'hnsw' is";
        " hyperparameter M.";
        "Note that some indices may not be able to return all the existing";
        "nearest neighbors" ],
      TA.Default (fun () -> Interfaiss.Type.to_string Defaults.neighbors_index_type),
      (fun _ ->
        Set_neighbors_index_type (TA.get_parameter () |> Interfaiss.Type.of_string) |> List.accum Parameters.program);
    [ "--neighbors-summarize-at-most"; "--neighbors-in-summary" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of closest sequences to be printed";
        "when summarizing nearest neighbors.";
        "Note that more might be printed anyway in case of ties.";
        "The statistics in the summary will be computed on all the neighbors";
        "explored according to the policy specified by option";
        "'--neighbors-guard-policy'" ],
      TA.Default (fun () -> KeepAtMost.to_string Defaults.neighbors_keep_at_most),
      (fun _ ->
        Set_neighbors_keep_at_most (TA.get_parameter () |> KeepAtMost.of_string) |> List.accum Parameters.program);
    [ "--neighbors-guard-policy"; "--neighbors-exploration-policy" ],
      Some "'times('<float_no_less_than_one>')'|'plus(<non-negative_integer>)'",
      [ "set the number of nearest neighbors to be explored";
        "when summarizing them.";
        "Note that this is greater than or equal to the number of neighbors";
        "specified with option '--neighbors-summarize-at-most'.";
        "Calling the latter n,";
        " policy 'times('m')'";
        "will explore m*n nearest neighbors, while";
        " policy 'plus('m')'";
        "will explore m+n nearest neighbors.";
        "The additional neighbors explored are not printed,";
        "but used to compute overall statistics" ],
      TA.Default (fun () -> Twisted.NeighborsPolicy.to_string Defaults.neighbors_guard_policy),
      (fun _ ->
        Set_neighbors_guard_policy
          (TA.get_parameter () |> Twisted.NeighborsPolicy.of_string) |> List.accum Parameters.program);
    [ "-n"; "--summarize-neighbors"; "--find-and-summarize-neighbors" ],
      Some "<twisted_binary_file_prefix> <summary_file_prefix>",
      [ "for each vector present in the twisted register, find nearest neighbors";
        "among the vectors present in the specified twisted binary file";
        " (which must have extension '.KPopTwisted' unless file is '/dev/*')";
        "using the current metric function, distance function and normalization;";
        "summarize distances and write the result to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.KPopSummary.txt' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ ->
        let twisted_prefix = TA.get_parameter () in
        Summary_from_twisted_neighbors (twisted_prefix, TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Experimental actions - They may be removed from future versions:" ];
    [ "--precision-for-splits" ],
      Some "<positive_integer>",
      [ "set how many precision digits should be used when outputting splits";
        "in plain-text format" ],
      TA.Default (fun () -> string_of_int Defaults.precision_splits),
      (fun _ -> Set_precision_splits (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "--splits-algorithm" ],
      Some "'gaps'|'centroids'",
      [ "algorithm to use when computing splits from embeddings" ],
      TA.Default (fun () -> Twisted.SplitsAlgorithm.to_string Defaults.splits_algorithm),
      (fun _ ->
        Set_splits_algorithm (TA.get_parameter () |> Twisted.SplitsAlgorithm.of_string)
        |> List.accum Parameters.program);
    [ "--splits-at-most"; "--splits-keep-at-most" ],
      Some "<positive_integer>|'all'",
      [ "set the maximum number of phylogenetic splits to be kept";
        "when generating them from embeddings" ],
      TA.Default (fun () -> string_of_int Defaults.splits_keep_at_most),
      (fun _ -> Set_splits_keep_at_most (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-S"; "--splits"; "--compute-splits"; "--twisted-to-splits" ],
      Some "<phylosplits_tabular_file_prefix>",
      [ "compute phylogenetic splits";
        "from the vectors present in the twisted register";
        "using the current metric function, distance function and normalization.";
        "The result will be written to the specified tabular file.";
        "File extension is automatically assigned";
        " (will be '.PhyloSplits' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Splits_from_twisted (TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately" ];
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (fun () -> "quiet execution"),
      (fun _ -> Parameters.verbose := true);
    [ "-V"; "--version" ],
      None,
      [ "print version and exit" ],
      TA.Optional,
      (fun _ -> Printf.printf "%s\n%!" info.version; exit 0);
    (* Hidden option to profile twisting *)
    [ "--debug-twisting" ], None, [], TA.Optional, (fun _ -> Parameters.debug_twisting := true);
    (* Hidden option to emit help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    (* Hidden option to print exception backtrace *)
    [ "-x"; "--print-exception-backtrace" ], None, [], TA.Optional, (fun _ -> Printexc.record_backtrace true);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 0)
  ];
  let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;
  if !Parameters.verbose then
    TA.header ();
  (* We perform a dry run of the program to detect possible errors *)
  let twister_loaded = ref false
  and distance = ref Defaults.distance and distance_normalize = ref Defaults.distance_normalize in
  List.iter
    (function
      | Empty _ ->
        ()
      | Binary_to_register (Twister, _) | Tables_to_register (Twister, _) ->
        twister_loaded := true
      | Binary_to_register (Twisted, _) | Tables_to_register (Twisted, _)
      | Add_binary_to_twisted _ ->
        ()
      | Twist_database _ ->
        (* A twister must have been loaded to twist spectra *)
        if not !twister_loaded then
          TA.parse_error
            "Option '-k'/'-s' requires a twister in the twister register!"
      | Register_to_binary _
      | Set_precision_tables _ | Set_precision_splits _ | Register_to_tables _
      | Set_metric _ ->
        ()
      | Set_distance dist ->
        distance := dist
      | Set_distance_normalize norm ->
        distance_normalize := norm
      | Embeddings_from_twisted _ | Splits_from_twisted _
      | Summary_from_twisted_binary _ | Summary_from_twisted_neighbors _ ->
        begin match !distance, !distance_normalize with
        | Space.Distance.Cosine, false | Angle, false ->
          TA.parse_error "Distances 'cosine' and 'angle' require embeddings to be normalized"
        | Cosine, true | Angle, true | Euclidean, _ | Minkowski _, _ ->
          ()
        end
      | Set_splits_algorithm _ | Set_splits_keep_at_most _ | Set_summary_keep_at_most _
      | Set_neighbors_keep_at_most _ | Set_neighbors_guard_policy _ | Set_neighbors_index_type _ ->
        ())
    program;
  (* These are the registers available to the program *)
  let twister = ref Twister.empty and twisted = ref Twisted.empty and metric = ref Defaults.metric
  and distance = ref Defaults.distance and distance_normalize = ref Defaults.distance_normalize
  and splits_keep_at_most = ref Defaults.splits_keep_at_most and splits_algorithm = ref Defaults.splits_algorithm
  and summary_keep_at_most = ref Defaults.summary_keep_at_most
  and neighbors_keep_at_most = ref Defaults.neighbors_keep_at_most
  and neighbors_guard_policy = ref Defaults.neighbors_guard_policy
  and neighbors_index_type = ref Defaults.neighbors_index_type
  and precision_tables = ref Defaults.precision_tables and precision_splits = ref Defaults.precision_splits in
  let twisted_of_binary = Twisted.of_binary ~verbose:!Parameters.verbose
  and twisted_of_files = Twisted.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose in
  try
    List.iter
      (function
        | Empty Twister ->
          twister := Twister.empty
        | Empty Twisted ->
          twisted := Twisted.empty
        | Binary_to_register (Twister, prefix) ->
          twister := Twister.of_binary ~verbose:!Parameters.verbose prefix
        | Binary_to_register (Twisted, prefix) ->
          twisted := twisted_of_binary prefix
        | Tables_to_register (Twister, prefix) ->
          twister := Twister.of_files ~threads:!Parameters.threads ~verbose:!Parameters.verbose prefix
        | Tables_to_register (Twisted, prefix) ->
          twisted := twisted_of_files prefix
        | Add_binary_to_twisted prefix ->
          twisted := twisted_of_binary prefix |> Twisted.merge_rowwise !twisted
        | Twist_database fname ->
          twisted :=
            Twister.add_twisted_from_database
              ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              ~debug:!Parameters.debug_twisting !twister !twisted fname
        | Register_to_binary (Twister, prefix) ->
          Twister.to_binary ~verbose:!Parameters.verbose !twister prefix
        | Register_to_binary (Twisted, prefix) ->
          Twisted.to_binary ~verbose:!Parameters.verbose !twisted prefix
        | Set_precision_tables prec ->
          precision_tables := prec
        | Set_precision_splits prec ->
          precision_splits := prec
        | Register_to_tables (Twister, prefix) ->
          Twister.to_files
            ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !twister prefix
        | Register_to_tables (Twisted, prefix) ->
          Twisted.to_files
            ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !twisted prefix
        | Set_metric metr ->
          metric := metr
        | Set_distance dist ->
          distance := dist
        | Set_distance_normalize norm ->
          distance_normalize := norm
        | Embeddings_from_twisted prefix ->
          let res =
            Twisted.to_embeddings
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance !metric !twisted in
          Matrix.to_file
            ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            res prefix
        | Set_splits_algorithm algo ->
          splits_algorithm := algo
        | Set_splits_keep_at_most kam ->
          splits_keep_at_most := kam
        | Splits_from_twisted prefix ->
          let res =
            Twisted.get_splits
              ~normalize:!distance_normalize ~threads:!Parameters.threads ~verbose:!Parameters.verbose
              !distance !metric !splits_algorithm !splits_keep_at_most !twisted in
          Trees.Splits.to_file ~precision:!precision_splits res prefix
        | Set_summary_keep_at_most kam ->
          summary_keep_at_most := kam
        | Summary_from_twisted_binary (prefix_in, prefix_out, output_distance_matrix) ->
          Twisted.summarize_distances_rowwise
            ~normalize:!distance_normalize ~keep_at_most:!summary_keep_at_most ~output_distance_matrix
            ~precision:!precision_tables ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !distance !metric (twisted_of_binary prefix_in) !twisted prefix_out
        | Set_neighbors_keep_at_most nhm ->
          neighbors_keep_at_most := nhm
        | Set_neighbors_guard_policy gp ->
          neighbors_guard_policy := gp
        | Set_neighbors_index_type it ->
          neighbors_index_type := it
        | Summary_from_twisted_neighbors (prefix_in, prefix_out) ->
          Twisted.summarize_neighbors
            ~normalize:!distance_normalize ~how_many:!neighbors_keep_at_most
            ~policy:!neighbors_guard_policy ~index_type:!neighbors_index_type
            ~threads:!Parameters.threads ~verbose:!Parameters.verbose
            !metric (twisted_of_binary prefix_in) !twisted prefix_out)
      program
  with
  | Exception.E (Exception.Kind.Initialize, _, _) | Exception.E (Exception.Kind.IO_Format, _, _) as e ->
    Exception.to_string e |> String.TermIO.red |> Printf.eprintf "(%s): FATAL: %s\n%!" __FUNCTION__
  | exc ->

(* TODO: WE SHOULD EXCLUDE THE CASE OF BROKEN PIPE *)

    Printf.peprintf "(%s): %s\n%!" __FUNCTION__
      ("FATAL: Uncaught exception: " ^ Printexc.to_string exc |> String.TermIO.red);
    Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
    Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__;
    Printexc.print_backtrace stderr

