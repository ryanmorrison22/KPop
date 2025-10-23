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
    exception Invalid_field of int * int
    (* Field number is 1-based! *)
    let extract n_field s =
      let numbers =
        try
          (* In principle this part should be safe *)
          Str.full_split regexp s
            |> List.filter_map
                 (function
                   | Str.Delim s -> Some (float_of_string s |> ceil |> int_of_float)
                   | Text _ -> None)
            |> Array.of_list
        with _ ->
          assert false in
      let n = Array.length numbers in
      if n_field > n then
        Invalid_field (n_field, n) |> raise
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

module Parameters =
  struct
    let option_l_or_L = ref false
    let content = KMI.Content.of_string "ds-DNA" |> ref
    let hasher = KMI.Hasher.K_mers 12 |> ref
    let weight_field = ref 0
    let max_results_size = ref 16777216 (* Or: 4^12 *)
    let inputs = ref []
    let label = ref ""
    let output = ref ""
    (*let threads = Tools.Parallel.get_nproc () |> ref*)
    let verbose = ref false
  end

let info = {
  Tools.Argv.name = "KPopCount";
  version = "21";
  date = "23-Oct-2025"
} and authors = [
  "2017-2025", "Paolo Ribeca", "paolo.ribeca@gmail.com"
]

let () =
  let module TA = Tools.Argv in
  TA.set_header (info, authors, [ BiOCamLib.Info.info; KPop.Info.info ]);
  TA.set_synopsis "-l <output_vector_label>|-L [OPTIONS]";
  TA.parse [
    TA.make_separator "Algorithmic parameters";
    [ "-k"; "--k-mer-size"; "--k-mer-length" ],
      Some "<positive_integer>",
      [ "set the hashing strategy to iteration over regular k-mers";
        "and specify the k-mer length to be used.";
        "Options '-k' and '-g' are mutually exclusive; if multiple are specified";
        "the last one will take effect" ],
      TA.Default
        (fun () ->
          match !Parameters.hasher with
          | K_mers _ -> KMI.Hasher.to_string !Parameters.hasher
          | Gapped _ -> "not used"),
      (fun _ -> Parameters.hasher := KMI.Hasher.K_mers (TA.get_parameter_int_pos ()));
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
          (fun () ->
            match !Parameters.hasher with
            | K_mers _ -> "not_used"
            | Gapped _ -> KMI.Hasher.to_string !Parameters.hasher),
      (fun _ ->
        let k = TA.get_parameter_int_pos () in
        let g = TA.get_parameter_int_pos () in
        Parameters.hasher := KMI.Hasher.Gapped (k, g));
    [ "--max-results-size" ],
      Some "<positive_integer>",
      [ "maximum number of k-mer hashes to be kept in memory at any given time.";
        "If more are present, the ones corresponding to the lowest cardinality";
        "will be removed from memory and printed out, and there will be";
        "repeated hashes in the output" ],
      TA.Default (fun () -> string_of_int !Parameters.max_results_size),
      (fun _ -> Parameters.max_results_size := TA.get_parameter_int_pos ());
    TA.make_separator "Input/Output";
    [ "-c"; "--content" ],
      Some "'ss-DNA'|'single-stranded-DNA'|'ds-DNA'|'double-stranded-DNA'|'protein'|FULL",
      [ "how file contents should be interpreted.";
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
      TA.Default (fun _ -> KMI.Content.to_string !Parameters.content),
      (fun _ -> Parameters.content := TA.get_parameter () |> KMI.Content.of_string);
    [ "-f"; "--fasta" ],
      Some "<fasta_file_name>",
      [ "FASTA input file containing sequences.";
        "You can specify more than one FASTA input, but not FASTA and FASTQ inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ -> Files.Type.FASTA (TA.get_parameter ()) |> List.accum Parameters.inputs);
    [ "-s"; "--single-end" ],
      Some "<fastq_file_name>",
      [ "FASTQ input file containing single-end sequencing reads";
        "You can specify more than one FASTQ input, but not FASTQ and FASTA inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ -> SingleEndFASTQ (TA.get_parameter ()) |> List.accum Parameters.inputs);
    [ "-p"; "--paired-end" ],
      Some "<fastq_file_name1> <fastq_file_name2>",
      [ "FASTQ input files containing paired-end sequencing reads";
        "You can specify more than one FASTQ input, but not FASTQ and FASTA inputs";
        "at the same time. Contents are expected to be homogeneous across inputs" ],
      TA.Optional,
      (fun _ ->
        let name1 = TA.get_parameter () in
        let name2 = TA.get_parameter () in
        PairedEndFASTQ (name1, name2) |> List.accum Parameters.inputs);
    [ "-w"; "--weights"; "--weights-from-sequence-names" ],
      Some "<non_negative_integer>",
      [ "given the index n specified as a parameter, extract the n-th number";
        "from each sequence name and weigh the corresponding sequence accordingly.";
        "Indices are 1-based; a value of 0 disables weighting.";
        "If the weight is a float number, the ceiling of such number will be used" ],
      TA.Default
        (fun _ ->
          if !Parameters.weight_field = 0 then
            "do not weigh"
          else
            string_of_int !Parameters.weight_field),
      (fun _ -> Parameters.weight_field := TA.get_parameter_int_non_neg ());
    [ "-l"; "--label" ],
      Some "<output_vector_label>",
      [ "label to be given to the k-mer spectrum in the output file.";
        "It must not contain double quote '\"' characters.";
        "Either option '-l' or option '-L' is mandatory" ],
      TA.Optional,
      (fun _ ->
        Parameters.option_l_or_L := true;
        Parameters.label :=
          let res = TA.get_parameter () in
          try
            Matrix.Base.strip_external_quotes_and_check res
          with Matrix.Base.Quotes_in_name _ ->
            TA.parse_error "Spectrum labels must not contain quotes";
            assert false); (* To keep the compiler happy *)
    [ "-L"; "--one-spectrum-per-sequence" ],
      None,
      [ "output one spectrum per input sequence, using the sequence name as label.";
        "Sequence names must not contain double quote '\"' characters.";
        "Either option '-l' or option '-L' is mandatory" ],
      TA.Optional,
      (fun _ -> Parameters.option_l_or_L := true);
    [ "-o"; "--output" ],
      Some "<output_file_prefix>",
      [ "prefix of the generated output file";
        " (will be given extension '.KPopSpectra.txt' unless file is '/dev/*')" ],
      TA.Default (fun () -> if !Parameters.output = "" then "<stdout>" else !Parameters.output),
      (fun _ -> Parameters.output := TA.get_parameter () |> KMerDB.Spectra.make_filename);
    TA.make_separator "Miscellaneous";
(*
    [ "-t"; "-T"; "--threads" ],
      Some "<computing_threads>",
      [ "number of concurrent computing threads to be spawned";
        " (default automatically detected from your configuration)" ],
      TA.Default (fun () -> string_of_int !Parameters.threads),
      (fun _ -> Parameters.threads := TA.get_parameter_int_pos ());
*)
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
    (* Hidden option to output help in markdown format *)
    [ "--markdown" ], None, [], TA.Optional, (fun _ -> TA.markdown (); exit 0);
    [ "-h"; "--help" ],
      None,
      [ "print syntax and exit" ],
      TA.Optional,
      (fun _ -> TA.usage (); exit 1)
  ];
  if not !Parameters.option_l_or_L then
    TA.parse_error "One of options '-l' and '-L' is mandatory";
  if !Parameters.verbose then
    TA.header ();
  Parameters.inputs := List.rev !Parameters.inputs;
  if !Parameters.inputs <> [] then begin
    let is_format_fasta = ref false and store = ref Files.ReadsIterate.empty in
    List.iteri
      (fun i input ->
        if i = 0 then
          is_format_fasta := begin
            match input with
            | Files.Type.FASTA _ -> true
            | SingleEndFASTQ _ | PairedEndFASTQ _ -> false
            | _ -> assert false
          end
        else
          if begin
            match input with
            | FASTA _ -> not !is_format_fasta
            | SingleEndFASTQ _ | PairedEndFASTQ _ -> !is_format_fasta
            | _ -> assert false
          end then
            TA.parse_error "You cannot process FASTA and FASTQ inputs together";
        store := Files.ReadsIterate.add_from_files !store [| input |])
      !Parameters.inputs;
    let output =
      if !Parameters.output = "" then
        stdout
      else
        open_out !Parameters.output in
    (* Header with label *)
    if !Parameters.label <> "" then
      Printf.fprintf output "\t%s\n" !Parameters.label;
    let reads_cntr = ref 0
    and k_mer_iterator =
      KMI.make ~max_results_size:!Parameters.max_results_size ~verbose:!Parameters.verbose
        !Parameters.content !Parameters.hasher (Printf.fprintf output "%s\t%d\n") in
    (* Note that linting is done automatically at a lower level by KMerIterator
        depending on the sequence type, so we disable it here *)
    Files.ReadsIterate.iter ~linter:Sequences.Lint.none ~verbose:false
      (fun _ segm_id read ->
        let read_tag = Matrix.Base.strip_external_quotes_and_check read.tag in
        (* If no global label has been specified, we output one per sequence *)
        if !Parameters.label = "" then
          Printf.fprintf output "\t%s\n" read_tag;
        let weight =
          if !Parameters.weight_field = 0 then
            1
          else
            CoverageFromName.extract !Parameters.weight_field read_tag in
        k_mer_iterator ~weight read.seq;
        if !Parameters.verbose && !reads_cntr mod 1_000 = 0 then
          Printf.eprintf "%s\r(%s): Added and hashed %d %s%!" String.TermIO.clear __FUNCTION__
            !reads_cntr (String.pluralize_int "read" !reads_cntr);
        if segm_id = 0 then
          incr reads_cntr)
      !store;
    if !Parameters.verbose then
      Printf.eprintf "%s\r(%s): Added and hashed %d %s.\n%!" String.TermIO.clear __FUNCTION__
        !reads_cntr (String.pluralize_int "read" !reads_cntr);
    (*Printf.eprintf "Times: (encode=%g, trie=%g, array=%g, accumulate=%g)\n%!"
      (Tools.Timer.read "KMers.Iterator.Encoder:encode")
      (Tools.Timer.read "KMers.Iterator.Encoder:trie")
      (Tools.Timer.read "KMers.Iterator.Encoder:array")
      (Tools.Timer.read "KMers.Iterator.Encoder:accumulate");*)
    close_out output
  end

