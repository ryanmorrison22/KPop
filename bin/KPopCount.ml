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

module CoverageFromName =
  struct
    let regexp = Str.regexp "[+-]?\\(\\([0-9]+[.]?[0-9]*\\)\\|\\([.][0-9]+\\)\\)\\([eE][-+]?[0-9]+\\)?"
    (* Field number is 1-based! *)
    let extract n_field s =
      let numbers =
        try
          (* In principle this part should be safe *)
          Str.full_split regexp s
            |> List.filter_map
                 (function
                   | Str.Delim s -> Some (float_of_string s)
                   | Text _ -> None)
            |> Array.of_list
        with _ ->
          assert false in
      let n = Array.length numbers in
      if n_field > n then
        Exception.raise __FUNCTION__ IO_Format (Printf.sprintf "Invalid field %d (found %d)" n_field n)
      else
        numbers.(n_field - 1)
  end

module KMerIterator =
  struct
    include KMers.Iterator
    module Hasher =
      struct
        include KMers.Iterator.Hasher
        (* We overwrite the stock library function to better suit this program *)
        let to_string = function
        | K_mers k -> Printf.sprintf "continuous k-mers of size %d" k
        | Gapped (k, g) -> Printf.sprintf "gapped k-mers of size %d (%d+%d+%d)" (2*k+g) k g k
      end
  end
module KMI = KMerIterator

type to_do_t =
  | Empty
  | Of_file of string
  | Set_content of KMI.Content.t
  | Set_hasher of KMI.Hasher.t
  | Set_weight_extractor of int
  | Set_max_results_size of int
  (* The string is the label - if empty, each sequence names is treated as a label *)
  | Add_sequences of Files.Type.t * string
  | To_file of string

module Defaults =
  struct
    let content = KMI.Content.of_string "ds-DNA"
    let hasher = KMI.Hasher.K_mers 12
    let weight_field = 0
    let max_results_size = 16777216 (* Or: 4^12 *)
    (*let threads = Tools.Parallel.get_nproc ()*)
    let verbose = false
  end

module Parameters =
  struct
    let program = ref []
    (*let threads = ref Defaults.threads*)
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPopCount";
  version = "23";
  date = "09-Nov-2025"
} and authors = [
  "2017-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "[ACTIONS]";
  TA.parse [
    TA.make_separator_multiline [ "Actions."; "They are executed delayed and in order of specification." ];
    TA.make_separator_multiline [ ""; "Input/Output of spectra databases:" ];
    [ "-0"; "--zero"; "--empty" ],
      None,
      [ "load an empty database into the register" ],
      TA.Optional,
      (fun _ -> Empty |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load into the register the database present in the specified file";
        " (which must have extension '.KPopCounter' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_file (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "save the database present in the register to the specified file";
        " (which will be given extension '.KPopCounter' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> To_file (TA.get_parameter ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Algorithmic parameters:" ];
    [ "-k"; "--k-mer-size"; "--k-mer-length" ],
      Some "<positive_integer>",
      [ "set the hashing strategy to iteration over regular k-mers";
        "and specify the k-mer length to be used.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> KMI.Hasher.to_string Defaults.hasher
          | Gapped _ -> "not used") |> Fun.const),
      (fun _ -> Set_hasher (KMI.Hasher.K_mers (TA.get_parameter_int_pos ())) |> List.accum Parameters.program);
    [ "-g"; "--gapped-k-mer-sizes"; "--gapped-k-mer-lengths" ],
      Some "BLOCK_SIZE GAP_SIZE",
      [ "where";
        " BLOCK-SIZE := <positive_integer>";
        " GAP-SIZE := <positive_integer>";
        "Set the hashing strategy to iteration over symmetrical gapped k-mers";
        "(having a BLOCK-GAP-BLOCK structure, with BLOCKs of the same size) and";
        "specify their geometry in terms of BLOCK and GAP sizes, respectively.";
        "For instance, option";
        " '-g 5 1'";
        "will iterate on all existing k-mers of size 11 (5+1+5) and not take the";
        "central nucleotide into account for the purpose of computing the hash.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified";
        "the last one will take effect" ],
      TA.Default
        ((match Defaults.hasher with
          | K_mers _ -> "not_used"
          | Gapped _ -> KMI.Hasher.to_string Defaults.hasher) |> Fun.const),
      (fun _ ->
        let k = TA.get_parameter_int_pos () in
        let g = TA.get_parameter_int_pos () in
        Set_hasher (KMI.Hasher.Gapped (k, g)) |> List.accum Parameters.program);
    [ "--max-results-size" ],
      Some "<positive_integer>",
      [ "maximum number of k-mer hashes to be kept in memory at any given time.";
        "If more are present, the ones corresponding to the lowest cardinality";
        "will be removed from memory and printed out, and there will be";
        "repeated hashes in the output" ],
      TA.Default (string_of_int Defaults.max_results_size |> Fun.const),
      (fun _ -> Set_max_results_size (TA.get_parameter_int_pos ()) |> List.accum Parameters.program);
    [ "-c"; "--content" ],
      Some "'ss-DNA'|'single-stranded-DNA'|'ds-DNA'|'double-stranded-DNA'|'protein'|FULL",
      [ "set how contents of following input files should be interpreted.";
        "When content is 'ss-DNA', 'protein' or 'text', only the sequence is hashed;";
        "when content is 'ds-DNA', both sequence and reverse complement are hashed.";
        "'ss-DNA' prevents automatic matching of reverse-complemented sequences;";
        "use it only when comparing a set of single, homogeneus sequences.";
        "These are shortcuts for the full form of this option, which is defined as";
        " FULL := 'DNA('STRANDEDNESS','CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION')'";
        "       | 'protein('UNKNOWN_CHAR_ACTION')'";
        "       | 'text('CASE_SENSITIVITY','UNKNOWN_CHAR_ACTION',";
        "               '<dictionary_file_name>')'";
        "where";
        " STRANDEDNESS := 'ss'|'single-stranded'|'ds'|'double-stranded'";
        " CASE_SENSITIVITY := 'ci'|'case-insensitive'|'cs'|'case-sensitive'";
        " UNKNOWN_CHAR_ACTION := 'split'|'ignore'|'error'";
        "If 'case-insensitive' is specified, DNA/protein sequences are converted to";
        "uppercase characters, while text sequences (and dictionary entries)";
        "are converted to lowercase characters.";
        "UNKNOWN_CHAR_ACTION decides what happens when an unknown character is found";
        "in the input. Option 'split' (the default for DNA/protein sequences) splits";
        "the input sequence and skips unknown characters (for instance N/X) whenever";
        "they are encountered; option 'ignore' (the default for text) silently skips";
        "unknown characters (for instance whitespace); option 'error' causes the";
        "program to abort.";
        "If a dictionary file is specified, each of its lines is interpreted";
        "as a different dictionary entry/token" ],
      TA.Default (KMI.Content.to_string Defaults.content |> Fun.const),
      (fun _ -> Set_content (TA.get_parameter () |> KMI.Content.of_string) |> List.accum Parameters.program);
    [ "-w"; "--weights"; "--weights-from-sequence-names" ],
      Some "<non_negative_integer>",
      [ "given the index n specified as a parameter, extract the n-th number";
        "from each sequence name and weigh the corresponding sequence accordingly.";
        "Indices are 1-based; a value of 0 disables weighting.";
        "If no such field exists, the program will fail.";
        "If the weight is a float number, the ceiling of such number will be used" ],
      TA.Default
        ((if Defaults.weight_field = 0 then
            "do not weigh"
          else
            string_of_int Defaults.weight_field) |> Fun.const),
      (fun _ -> Set_weight_extractor (TA.get_parameter_int_non_neg ()) |> List.accum Parameters.program);
    TA.make_separator_multiline [ ""; "Input/Output of sequences for processing:" ];
    [ "-f"; "--fasta" ],
      Some "<fasta_file_name> <label>",
      [ "process the sequences contained in the specified FASTA input file.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let path = TA.get_parameter () in
        let label = TA.get_parameter () in
        Add_sequences (Files.Type.FASTA path, label) |> List.accum Parameters.program);
    [ "-s"; "--single-end" ],
      Some "<fastq_file_name> <label>",
      [ "process the sequences contained in the specified FASTQ input file";
        "containing single-end sequencing reads.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let path = TA.get_parameter () in
        let label = TA.get_parameter () in
        Add_sequences (SingleEndFASTQ path, label) |> List.accum Parameters.program);
    [ "-p"; "--paired-end" ],
      Some "<fastq_file_name1> <fastq_file_name2> <label>",
      [ "process the sequences contained in the specified FASTQ input file";
        "containing paired-end sequencing reads.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let path1 = TA.get_parameter () in
        let path2 = TA.get_parameter () in
        let label = TA.get_parameter () in
        Add_sequences (PairedEndFASTQ (path1, path2), label) |> List.accum Parameters.program);
    [ "-t"; "--tabular" ],
      Some "<fasta_file_name> <label>",
      [ "process the sequences contained in the specified tabular input file.";
        "If a label is specified, the hashes extracted from all the sequences";
        "are collected into one spectrum having the label as name; if the label is";
        "empty, each sequence is turned into one separate spectrum and the";
        "sequence name is used as label. Label and sequence names must not contain";
        "double quote '\"' characters.";
        "While you can specify several inputs possibly having different formats,";
        "contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let path = TA.get_parameter () in
        let label = TA.get_parameter () in
        Add_sequences (Files.Type.Tabular path, label) |> List.accum Parameters.program);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately" ];
(*
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (string_of_int Defaults.threads |> Fun.const),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
*)
    [ "-v"; "--verbose" ],
      None,
      [ "set verbose execution" ],
      TA.Default (Fun.const "quiet execution"),
      (fun _ -> Parameters.verbose := true);
    [ "-V"; "--version" ],
      None,
      [ "print version and exit" ],
      TA.Optional,
      (fun _ -> Printf.printf "%s\n%!" info.version; exit 0);
    (* Hidden option to output help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  let program = List.rev !Parameters.program in
  if program = [] then begin
    TA.usage ();
    exit 0
  end;
  if !Parameters.verbose then
    TA.header ();
  (* These are the registers available to the program *)
  let db = KMerDB.make_empty () |> ref and content = ref Defaults.content and hasher = ref Defaults.hasher
  and weight_field = ref Defaults.weight_field and max_results_size = ref Defaults.max_results_size
  and reads_cntr = ref 0
  and catch_unexpected_end_of_output_file f =
    try
      f ()
    with End_of_file ->
      Exception.raise_unexpected_end_of_output __FUNCTION__ in
  try
    List.iter
      (function
        | Empty ->
          db := KMerDB.make_empty ()
        | Of_file prefix ->
          db := KMerDB.of_binary ~verbose:!Parameters.verbose prefix
        | Set_content c ->
          content := c
        | Set_hasher h ->
          hasher := h
        | Set_weight_extractor wf ->
          weight_field := wf
        | Set_max_results_size mrs ->
          max_results_size := mrs
        | Add_sequences (input, label) ->
          let current = ref label in
          let k_mer_iterator =
            KMI.make ~max_results_size:!max_results_size ~verbose:!Parameters.verbose
              !content !hasher (fun hash n -> db := KMerDB.add_counts !db !current hash n) in
          (* Note that linting is done automatically at a lower level by KMerIterator
              depending on the sequence type, so we disable it here *)
          Files.ReadsIterate.iter ~linter:Sequences.Lint.none ~verbose:false
            (fun _ segm_id read ->
              (* If no global label has been specified, we output one per sequence *)
              if label = "" then
                current := read.tag;
              let weight =
                if !weight_field = 0 then
                  1.
                else
                  CoverageFromName.extract !weight_field read.tag in
              k_mer_iterator ~weight read.seq;
              if !Parameters.verbose && !reads_cntr mod 1_000 = 0 then
                Printf.eprintf "%s\r(%s): Added and hashed %d %s%!" String.TermIO.clear __FUNCTION__
                  !reads_cntr (String.pluralize_int "read" !reads_cntr);
              if segm_id = 0 then
                incr reads_cntr)
            (Files.ReadsIterate.add_from_files Files.ReadsIterate.empty [| input |]);
          if !Parameters.verbose then
            Printf.eprintf "%s\r(%s): Added and hashed %d %s.\n%!" String.TermIO.clear __FUNCTION__
              !reads_cntr (String.pluralize_int "read" !reads_cntr);
        | To_file prefix ->
          catch_unexpected_end_of_output_file
            (fun () -> KMerDB.to_binary ~verbose:!Parameters.verbose !db prefix))
      program
    (*;Printf.eprintf "Times: (encode=%g, trie=%g, array=%g, accumulate=%g)\n%!"
      (Tools.Timer.read "KMers.Iterator.Encoder:encode")
      (Tools.Timer.read "KMers.Iterator.Encoder:trie")
      (Tools.Timer.read "KMers.Iterator.Encoder:array")
      (Tools.Timer.read "KMers.Iterator.Encoder:accumulate");*)
  with
  | Exception.E (Exception.Kind.Initialize, _, _) | Exception.E (Exception.Kind.IO_Format, _, _) as e ->
    TA.usage ();
    Exception.to_string e |> String.TermIO.red |> Printf.eprintf "(%s): FATAL: %s\n%!" __FUNCTION__
  | exc ->
    Printf.peprintf "(%s): %s\n%!" __FUNCTION__
      ("FATAL: Uncaught exception: " ^ Printexc.to_string exc |> String.TermIO.red);
    Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
    Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__;
    Printexc.print_backtrace stderr

