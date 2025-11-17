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

module Action =
  struct
    type t =
      | Empty
      | Of_file of string
      | Set_content of KMI.Content.t
      | Set_hasher of KMI.Hasher.t
      | Set_weight_extractor of int
      (* The string is the label - if empty, each sequence names is treated as a label *)
      | Add_sequences of (Files.Type.t * string) list
      | To_file of string
    let compact_add_sequences program file_type s =
      program :=
        match !program with
        | Add_sequences l :: tl ->
          Add_sequences ((file_type, s) :: l) :: tl
        | l ->
          Add_sequences [file_type, s] :: l
  end

module Defaults =
  struct
    let content = KMI.Content.of_string "ds-DNA"
    let hasher = KMI.Hasher.K_mers 12
    let weight_field = 0
    let threads = Processes.Parallel.get_nproc ()
    let verbose = false
  end

module Parameters =
  struct
    let program = ref []
    let threads = ref Defaults.threads
    let verbose = ref Defaults.verbose
  end

let info = {
  Tools.Argv.name = "KPopCount";
  version = "26";
  date = "16-Nov-2025"
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
      (fun _ -> Action.Empty |> List.accum Parameters.program);
    [ "-i"; "--input" ],
      Some "<binary_file_prefix>",
      [ "load into the register the database present in the specified file";
        " (which must have extension '.KPopSpectra' unless file is '/dev/*')" ],
      TA.Optional,
      (fun _ -> Of_file (TA.get_parameter ()) |> List.accum Parameters.program);
    [ "-o"; "--output" ],
      Some "<binary_file_prefix>",
      [ "save the database present in the register to the specified file";
        " (which will be given extension '.KPopSpectra' unless file is '/dev/*')" ],
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
        Action.compact_add_sequences Parameters.program (Files.Type.FASTA path) label);
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
        Action.compact_add_sequences Parameters.program (SingleEndFASTQ path) label);
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
        Action.compact_add_sequences Parameters.program (PairedEndFASTQ (path1, path2)) label);
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
        Action.compact_add_sequences Parameters.program (Files.Type.Tabular path) label);
    TA.make_separator_multiline [ "Miscellaneous options."; "They are set immediately" ];
    [ "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (string_of_int Defaults.threads |> Fun.const),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
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
  and weight_field = ref Defaults.weight_field in
  try
    List.iter
      (function
        | Action.Empty ->
          db := KMerDB.make_empty ()
        | Of_file prefix ->
          db := KMerDB.of_binary ~verbose:!Parameters.verbose prefix
        | Set_content c ->
          content := c
        | Set_hasher h ->
          hasher := h
        | Set_weight_extractor wf ->
          weight_field := wf
        | Add_sequences files ->
          (* Files will have been inverted during input compaction *)
          let files = Array.of_rlist files
          and file_cntr = ref 0 and spectrum_cntr = ref 0 in
          if !Parameters.threads = 1 || Array.length files = 1 then begin
            let current = -1 |> ref in
            let k_mer_iterator, k_mer_finalizer =
              KMI.make ~verbose:!Parameters.verbose !content !hasher
                (fun hash n -> KMerDB.add_counts_unsafe db !current hash n) in
            Array.iter
              (fun (input, label) ->
                if !Parameters.verbose then
                  Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s%!" String.TermIO.clear __FUNCTION__
                    !spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" !spectrum_cntr)
                    !file_cntr (String.pluralize_int "file" !file_cntr);
                incr file_cntr;
                (* Note that linting is done automatically at a lower level by KMerIterator
                    depending on the sequence type, so we disable it here *)
                Files.ReadsIterate.iter ~linter:Sequences.Lint.none ~verbose:false
                  (fun _ segm_id read ->
                    let weight =
                      if !weight_field = 0 then
                        1.
                      else
                        CoverageFromName.extract !weight_field read.tag in
                    k_mer_iterator ~weight read.seq;
                    (* If no global label has been specified, we output one spectrum per sequence *)
                    if label = "" then begin
                      current := KMerDB.add_empty_column_if_needed db read.tag;
                      k_mer_finalizer ();
                      if segm_id = 0 then
                        incr spectrum_cntr
                    end)
                  (Files.ReadsIterate.add_from_files Files.ReadsIterate.empty [| input |]);
                (* If a global label has been specified, we output one spectrum for the whole file *)
                if label <> "" then begin
                  current := KMerDB.add_empty_column_if_needed db label;
                  k_mer_finalizer ();
                  incr spectrum_cntr
                end)
              files
          end else begin
            let num_files = Array.length files and file_idx = ref 0
            and hashes = Tools.StackArray.create () |> ref and res = Tools.StackArray.create () in
            let k_mer_iterator, k_mer_finalizer =
              KMI.make ~verbose:!Parameters.verbose !content !hasher
                (fun hash n -> Tools.StackArray.push !hashes (hash, n)) in
            Processes.Parallel.process_stream_chunkwise
              (fun () ->
                assert (!file_idx <= num_files);
                if !file_idx = num_files then
                  raise End_of_file;
                let input, label = files.(!file_idx) in
                incr file_idx;
                Files.ReadsIterate.add_from_files Files.ReadsIterate.empty [| input |], label)
              (fun (iterator, label) ->
                Tools.StackArray.clear res;
                hashes := Tools.StackArray.create ();
                (* Note that linting is done automatically at a lower level by KMerIterator
                    depending on the sequence type, so we disable it here *)
                Files.ReadsIterate.iter ~linter:Sequences.Lint.none ~verbose:false
                  (fun _ segm_id read ->
                    let weight =
                      if !weight_field = 0 then
                        1.
                      else
                        CoverageFromName.extract !weight_field read.tag in
                    k_mer_iterator ~weight read.seq;
                    (* If no global label has been specified, we output one spectrum per sequence *)
                    if label = "" then begin
                      k_mer_finalizer ();
                      Tools.StackArray.push res (segm_id, read.tag, !hashes);
                      hashes := Tools.StackArray.create ()
                    end)
                  iterator;
                (* If a global label has been specified, we output one spectrum for the whole file *)
                if label <> "" then begin
                  k_mer_finalizer ();
                  Tools.StackArray.push res (0, label, !hashes)
                end;
                res)
              (fun res ->
                if !Parameters.verbose then
                  Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s%!" String.TermIO.clear __FUNCTION__
                    !spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" !spectrum_cntr)
                    !file_cntr (String.pluralize_int "file" !file_cntr);
                incr file_cntr;
                Tools.StackArray.riter
                  (fun (segm_id, label, hashes) ->
                    let current = KMerDB.add_empty_column_if_needed db label in
                    Tools.StackArray.riter
                      (fun (hash, n) ->
                        KMerDB.add_counts_unsafe db current hash n)
                      hashes;
                    if segm_id = 0 then
                      incr spectrum_cntr)
                  res)
              !Parameters.threads;
          end;
          if !Parameters.verbose then
            Printf.eprintf "%s\r(%s): Hashed and added %d %s from %d %s.\n%!" String.TermIO.clear __FUNCTION__
              !spectrum_cntr (String.pluralize_int ~plural:"spectra" "spectrum" !spectrum_cntr)
              !file_cntr (String.pluralize_int "file" !file_cntr)
        | To_file prefix ->
          Exception.catch_unexpected_end_of_output __FUNCTION__
            (fun () -> KMerDB.to_binary ~verbose:!Parameters.verbose !db prefix))
      program;
    if !Parameters.verbose && !Parameters.threads = 1 then
      Printf.eprintf "(%s): Timers: (iterate=%g (lint=%g, encode=%g, accumulate=%g), finalize=%g)\n%!" __FUNCTION__
        (Tools.Timer.read "KMers.Iterator.Iterator") (Tools.Timer.read "KMers.Iterator.Linter")
        (Tools.Timer.read "KMers.Iterator.Encoder") (Tools.Timer.read "KMers.Iterator.Accumulator")
        (Tools.Timer.read "KMers.Iterator.Finalizer")
  with e ->
    Exception.handle __FUNCTION__ TA.usage (fun () ->
      Printf.peprintf "(%s): This should not have happened - please contact <paolo.ribeca@gmail.com>\n%!" __FUNCTION__;
      Printf.peprintf "(%s): You might also wish to rerun me with option -x to get a full backtrace.\n%!" __FUNCTION__
    ) e

