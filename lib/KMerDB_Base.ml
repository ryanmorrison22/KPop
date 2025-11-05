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

(* We cannot open BiOCamLib here due to the ambiguity between BiOCamLib.Matrix and KPop.Matrix *)
module Numbers = BiOCamLib.Numbers
module Processes = BiOCamLib.Processes
module Tools = BiOCamLib.Tools
open BiOCamLib.Better

include (
  struct
    (* Utility function to resize arrays *)
    let _resize_t_array_ ?(is_buffer = true) length resize create_null nx ny a =
      let lx = Array.length a in
      let eff_nx =
        if is_buffer then begin
          if lx < nx then
            max nx (lx * 14 / 10)
          else
            lx
        end else
          nx
      and eff_ny =
        if is_buffer then begin
          if lx > 0 then begin
            (* We assume all bigarrays to have the same size *)
            let ly = length a.(0) in
            if ly < ny then
              max ny (ly * 14 / 10)
            else
              ly
          end else
            ny
        end else
          ny in
      (*Printf.eprintf "(%s): Resizing to (%d,%d) - asked (%d,%d)...\n%!" __FUNCTION__ eff_nx eff_ny nx ny;*)
      if eff_nx > lx then
        Array.append
          (* We need to provide the is_buffer argument like this because of type resolution *)
          (Array.map (resize ?is_buffer:(Some false) eff_ny) a)
          (Array.init (eff_nx - lx) (fun _ -> create_null eff_ny))
      else if eff_nx < lx then
        Array.map (resize ?is_buffer:(Some false) eff_ny) (Array.sub a 0 eff_nx)
      else (* eff_nx = lx *)
        if lx > 0 && eff_ny = length a.(0) then
          a
        else
          Array.map (resize ?is_buffer:(Some false) eff_ny) a
    module String =
      struct
        include String
        let resize_array ?(is_buffer = true) n a =
          Array.resize ~is_buffer n "" a
        (* *)
        let resize_array_array ?(is_buffer = true) =
          _resize_t_array_ ~is_buffer Array.length resize_array (fun l -> Array.make l "")
      end
    (* Counts are represented as 32-bit floats, with a couple additional modifications *)
    module CountBAVector =
      struct
        include Numbers.F32BAVector
        let resize ?(is_buffer = true) n =
          resize ~is_buffer ~fill_with:N.zero n
        let resize_array ?(is_buffer = true) =
          _resize_t_array_ ~is_buffer length resize (fun l -> make l N.zero)
      end
    (* *)
    type marshalled_t = {
      n_cols: int;
      n_rows: int;
      n_meta: int;
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      idx_to_meta_names: string array;
      meta: string array array;
      storage: CountBAVector.t array
    }
    and t = {
      core: marshalled_t;
      col_names_to_idx: int StringHashtbl.t;
      row_names_to_idx: int StringHashtbl.t;
      meta_names_to_idx: int StringHashtbl.t
    }
    (* *)
    let make_empty () = {
      core = {
        n_cols = 0;
        n_rows = 0;
        n_meta = 0;
        idx_to_col_names = [||];
        idx_to_row_names = [||];
        idx_to_meta_names = [||];
        meta = [||];
        storage = [||]
      };
      col_names_to_idx = StringHashtbl.create 16;
      row_names_to_idx = StringHashtbl.create 16;
      meta_names_to_idx = StringHashtbl.create 16
    }
    let output_summary ?(verbose = false) db =
      Printf.eprintf "[Spectrum labels (%d)]:" db.core.n_cols;
      Array.iteri
        (fun i s ->
          if i < db.core.n_cols then
            Printf.eprintf " '%s'" s)
        db.core.idx_to_col_names;
      Printf.eprintf "\n%!";
      if verbose then begin
        Printf.eprintf "[K-mer hashes (%d)]:" db.core.n_rows;
        Array.iteri
          (fun i s ->
            if i < db.core.n_rows then
              Printf.eprintf " '%s'" s)
          db.core.idx_to_row_names;
        Printf.eprintf "\n%!"
      end;
      Printf.eprintf "[Meta-data fields (%d)]:" db.core.n_meta;
      Array.iteri
        (fun i s ->
          if i < db.core.n_meta then
            Printf.eprintf " '%s'" s)
        db.core.idx_to_meta_names;
      Printf.eprintf "\n%!"
    (*  Utility functions  *)
    let invert_table a =
      let res = StringHashtbl.create (Array.length a) in
      Array.iteri (fun i name -> StringHashtbl.add res name i) a;
      res
    (* Makes space for a new sample named label *)
    let add_empty_column_if_needed db label =
      let n_cols = !db.core.n_cols in
      let aug_n_cols = n_cols + 1 in
      if StringHashtbl.mem !db.col_names_to_idx label |> not then begin
        StringHashtbl.add !db.col_names_to_idx label n_cols; (* THIS ONE CHANGES !db *)
        db := {
          !db with
          core = {
            !db.core with
            n_cols = aug_n_cols;
            (* We have to resize all the relevant containers *)
            idx_to_col_names = Array.append !db.core.idx_to_col_names [| label |];
            meta = String.resize_array_array ~is_buffer:true aug_n_cols !db.core.n_meta !db.core.meta;
            storage = CountBAVector.resize_array ~is_buffer:true aug_n_cols !db.core.n_rows !db.core.storage
          }
        }
      end
    (* Makes space for a new k-mer named hash *)
    let add_empty_row_if_needed db hash =
      if StringHashtbl.mem !db.row_names_to_idx hash |> not then begin
        let n_rows = !db.core.n_rows in
        let aug_n_rows = n_rows + 1 in
        StringHashtbl.add !db.row_names_to_idx hash n_rows; (* THIS ONE CHANGES !db *)
        db := {
          !db with
          core = {
            !db.core with
            n_rows = aug_n_rows;
            (* We have to resize all the relevant containers *)
            idx_to_row_names = begin
              let res = String.resize_array ~is_buffer:true aug_n_rows !db.core.idx_to_row_names in
              res.(n_rows) <- hash;
              res
            end;
            storage = CountBAVector.resize_array ~is_buffer:true !db.core.n_cols aug_n_rows !db.core.storage
          }
        }
      end
    (* *)
    exception Wrong_number_of_columns of int * int * int
    let add_metadata ?(verbose = false) db fname =
      let input = open_in fname and line_num = ref 0 in
      let header =
        input_line input
          |> String.Split.on_char_as_array '\t' |> Array.map Matrix.Base.strip_external_quotes_and_check in
      incr line_num;
      (* We add the names *)
      let missing = ref [] in
      Array.iteri
        (fun i name ->
          if i > 0 && StringHashtbl.mem db.meta_names_to_idx name |> not then
            List.accum missing name)
        header;
      let missing = Array.of_rlist !missing in
      let db = ref db
      and missing_len = Array.length missing in
      if missing_len > 0 then begin
        Array.iteri
          (fun i name -> !db.core.n_meta + i |> StringHashtbl.add !db.meta_names_to_idx name)
          missing;
        let n_meta = !db.core.n_meta + missing_len in
        db := {
          !db with
          core = {
            !db.core with
            n_meta;
            (* We have to resize all the relevant containers *)
            idx_to_meta_names = Array.append !db.core.idx_to_meta_names missing;
            meta = String.resize_array_array ~is_buffer:true !db.core.n_cols n_meta !db.core.meta
          }
        }
      end;
      let num_header_fields = Array.length header
      and meta_indices =
        Array.mapi
          (fun i name ->
            if i = 0 then
              -1
            else
              StringHashtbl.find !db.meta_names_to_idx name)
          header in
      begin try
        while true do
          let line =
            input_line input
              |> String.Split.on_char_as_array '\t' |> Array.map Matrix.Base.strip_external_quotes_and_check in
          incr line_num;
          (* A regular line. The first element is the spectrum name, the others the values of meta-data fields *)
          let l = Array.length line in
          if l <> num_header_fields then
            Wrong_number_of_columns (!line_num, l, num_header_fields) |> raise;
          add_empty_column_if_needed db line.(0);
          let col_idx = StringHashtbl.find !db.col_names_to_idx line.(0) in
          Array.iteri
            (fun i name_idx ->
              if i > 0 then
                !db.core.meta.(col_idx).(name_idx) <- line.(i))
            meta_indices;
          if verbose && !line_num mod 10 = 0 then
            Printf.eprintf "%s\r(%s): File '%s': Read %d lines%!"
              String.TermIO.clear __FUNCTION__ fname !line_num
        done
      with End_of_file ->
        close_in input;
        if verbose then
          Printf.eprintf "%s\r(%s): File '%s': Read %d lines\n%!"
            String.TermIO.clear __FUNCTION__ fname !line_num
      end;
      !db
    let add_counts db label hash counts =
      let db = ref db in
      add_empty_column_if_needed db label;
      add_empty_row_if_needed db hash;
      let col_idx = StringHashtbl.find !db.col_names_to_idx label
      and row_idx = StringHashtbl.find !db.row_names_to_idx hash in
      CountBAVector.(!db.core.storage.(col_idx).+(row_idx) <- N.of_int counts);
      !db
    let merge ?(verbose = false) db1 db2 =

      ignore (verbose, db1, db2);

      (*let area_1 = db_1.core.n_cols*)


      (* It's OK to merge DBs with different k-mer sets, but they must have identical metadata rows *)


      assert false

    (* *)
    let selected_from_regexps ?(verbose = false) db regexps =
      (* We iterate over the columns *)
      if verbose then
        Printf.eprintf "(%s): Selecting columns... [%!" __FUNCTION__;
      List.iter
        (fun (what, _) ->
          if verbose && what <> "" && StringHashtbl.find_opt db.meta_names_to_idx what = None then
            Printf.eprintf " (WARNING: Metadata field '%s' not found, no column will match)%!" what)
        regexps;
      let res = ref StringSet.empty in
      Array.iteri
        (fun n_col col_name ->
          if begin
            List.fold_left
              (fun start (what, regexp) ->
                start &&
                if what = "" then
                  (* Case of the label *)
                  Str.string_match regexp col_name 0
                else
                  match StringHashtbl.find_opt db.meta_names_to_idx what with
                  | None ->
                    false
                  | Some found ->
                    assert (db.core.idx_to_meta_names.(found) = what);
                    Str.string_match regexp db.core.meta.(n_col).(found) 0)
              true regexps
          end then
            res := StringSet.add col_name !res)
        db.core.idx_to_col_names;
      if verbose then
        StringSet.iter (Printf.eprintf " '%s'%!") !res;
      if verbose then
        Printf.eprintf " ] done.\n%!";
      !res
    let selected_negate db selection =
      StringSet.diff (Array.to_list db.core.idx_to_col_names |> StringSet.of_list) selection
    let remove_selected db selected =
      (* First, we compute the indices of the columns to be kept.
         We keep the same column order as in the original matrix *)
      let idxs = ref IntSet.empty in
      Array.iteri
        (fun col_idx col_name ->
          if StringSet.mem col_name selected |> not then
            idxs := IntSet.add col_idx !idxs)
        db.core.idx_to_col_names;
      let idxs = IntSet.elements_array !idxs in
      let n = Array.length idxs in
      let filter_array a = Array.init n (fun i -> a.(idxs.(i))) in
      let core =
        { db.core with
          n_cols = n;
          idx_to_col_names = filter_array db.core.idx_to_col_names;
          meta = filter_array db.core.meta;
          storage = filter_array db.core.storage } in
      { core = core;
        col_names_to_idx = invert_table core.idx_to_col_names;
        row_names_to_idx = invert_table core.idx_to_row_names;
        meta_names_to_idx = invert_table core.idx_to_meta_names }
    (* *)
    let archive_version = "2025-10-20"
    (* *)
    let make_filename_binary = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter"
    let to_binary ?(verbose = false) db prefix =
      let fname = make_filename_binary prefix in
      let output = open_out fname in
      if verbose then
        Printf.eprintf "(%s): Outputting DB to file '%s'...%!" __FUNCTION__ fname;
      output_value output "KPopCounter";
      output_value output archive_version;
      output_value output {
        db.core with
        (* We have to truncate all the containers *)
        idx_to_col_names = String.resize_array ~is_buffer:false db.core.n_cols db.core.idx_to_col_names;
        idx_to_row_names = String.resize_array ~is_buffer:false db.core.n_rows db.core.idx_to_row_names;
        idx_to_meta_names = String.resize_array ~is_buffer:false db.core.n_meta db.core.idx_to_meta_names;
        meta = String.resize_array_array ~is_buffer:false db.core.n_cols db.core.n_meta db.core.meta;
        storage = CountBAVector.resize_array ~is_buffer:false db.core.n_cols db.core.n_rows db.core.storage
      };
      close_out output;
      if verbose then
        Printf.eprintf " done.\n%!"
    let of_binary ?(verbose = false) prefix =
      let fname = make_filename_binary prefix in
      let input = open_in fname in
      if verbose then
        Printf.eprintf "(%s): Reading DB from file '%s'...%!" __FUNCTION__ fname;
      let which = (input_value input: string) in
      let version = (input_value input: string) in
      if which <> "KPopCounter" then
        Exception.raise __FUNCTION__ IO_Format
          (Printf.sprintf "Unexpected archive type (found '%s', expected 'KPopCounter')" which);
      if version <> archive_version then
        (* We are kind of misusing this function here *)
        Exception.raise_incompatible_archive_version __FUNCTION__ version archive_version;
      let core = (input_value input: marshalled_t) in
      close_in input;
      if verbose then
        Printf.eprintf " done.\n%!";
      { core = core;
        col_names_to_idx = invert_table core.idx_to_col_names;
        row_names_to_idx = invert_table core.idx_to_row_names;
        meta_names_to_idx = invert_table core.idx_to_meta_names }
  end: sig
    module CountBAVector:
      sig
        include module type of Numbers.F32BAVector
        val resize: ?is_buffer:bool -> int -> t -> t
      end
    (* Conceptually, each k-mer spectrum is stored as a column, even though in practice we store the transposed matrix -
        i.e., storage is a vector of spectra *)
    type marshalled_t = {
      n_cols: int; (* The number of spectra *)
      n_rows: int; (* The number of k-mers *)
      n_meta: int; (* The number of metadata fields *)
      (* We number rows, columns and metadata fields starting from 0 *)
      idx_to_col_names: string array;
      idx_to_row_names: string array;
      idx_to_meta_names: string array;
      (* *)
      meta: string array array; (* Dims = n_cols * n_meta *)
      storage: CountBAVector.t array (* Frequencies are stored as integers. Dims = n_cols * n_rows *)
    }
    and t = {
      core: marshalled_t;
      (* Inverted hashes for parsing *)
      col_names_to_idx: int StringHashtbl.t; (* Labels *)
      row_names_to_idx: int StringHashtbl.t; (* Hashes *)
      meta_names_to_idx: int StringHashtbl.t (* Metadata fields *)
    }
    val make_empty: unit -> t
    val merge: ?verbose:bool -> t -> t -> t
    val add_empty_column_if_needed: t ref -> string -> unit
    val add_empty_row_if_needed: t ref -> string -> unit
    (* Increase the counts for the specified sample and k-mer by the given amount *)
    val add_counts: t -> string -> string -> int -> t
    (* Add metadata - the first field must be the label *)
    exception Wrong_number_of_columns of int * int * int
    val add_metadata: ?verbose:bool -> t -> string -> t
    (* Select column names identified by regexps on metadata fields *)
    val selected_from_regexps: ?verbose:bool -> t -> (string * Str.regexp) list -> StringSet.t
    val selected_negate: t -> StringSet.t -> StringSet.t
    (* Remove spectra with the given labels *)
    val remove_selected: t -> StringSet.t -> t
    (* Output information about the contents *)
    val output_summary: ?verbose:bool -> t -> unit
    (* *)
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t (* Can fail due to archive version *)
  end
)

