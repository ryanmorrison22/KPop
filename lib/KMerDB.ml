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
    let ( .@() ) = CountBAVector.( .@() )
    let ( .@()<- ) = CountBAVector.( .@()<- )
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
    module ColOrRow =
      struct
        type t =
        | Col
        | Row
        let to_string = function
        | Col -> "column"
        | Row -> "row"
      end
    module Transformation =
      struct
        type statistics_t = {
          non_zero: int;
          min: int;
          max: int;
          sum: float;
          sum_log: float
        }
        type statistics_table_t = {
          col_stats: statistics_t array;
          row_stats: statistics_t array
        }
        type t =
          | Threshold of float
          | Power of float
          | CLR
          | Pseudo of float * bool
        let raise __FUNCTION__ kind s =
          Exception.raise __FUNCTION__ IO_Format (Printf.sprintf "%s transformation '%s'" kind s)
        let of_string_re = Str.regexp "[(,)]"
        let of_string s =
          let raise kind = raise __FUNCTION__ kind s in
          try
            match Str.full_split of_string_re s with
            | [ Text "threshold"; Delim "("; Text threshold; Delim ")" ] ->
              Threshold (float_of_string threshold)
            | [ Text "pow"; Delim "("; Text power; Delim ")" ]
            | [ Text "power"; Delim "("; Text power; Delim ")" ] ->
              Power (float_of_string power)
            | [ Text "binary" ] ->
              Power 0.
            | [ Text "clr" ] | [ Text "CLR" ] ->
              CLR
            | [ Text "pseudo"; Delim "("; Text power; Delim ","; Text quantize; Delim ")" ]
            | [ Text "pseudocounts"; Delim "("; Text power; Delim ","; Text quantize; Delim ")" ] ->
              let power = float_of_string power in
              let res = Pseudo (power, bool_of_string quantize) in
              if power < 0. then
                raise "Invalid";
              res
            | _ ->
              raise "Unknown"
          with _ ->
            raise "Unknown"
        let to_string = function
          | Threshold threshold ->
            Printf.sprintf "threshold(%g)" threshold
          | Power power ->
            Printf.sprintf "power(%g)" power
          | CLR ->
            "clr"
          | Pseudo (power, quantize) ->
            Printf.sprintf "pseudocounts(%g,%b)" power quantize
        let [@warning "-27"] compute ~which ~col_num ~col_stats ~row_num ~row_stats counts =
          match which with
          | Threshold threshold ->
            let threshold =
              if threshold < 1. then
                threshold *. col_stats.sum
              else
                threshold in
            if counts >= threshold then
              counts
            else
              0.
          | Power 1. -> (* Optimisation *)
            counts
          | Power power ->
            counts ** power
          | CLR ->
            (float_of_int col_stats.max +. 1.) *. log (counts +. 1.)
          | Pseudo (power, quantize) ->
            if power < 0. then
              to_string which |> raise __FUNCTION__ "Invalid";
            let v =
              if power = 0. then
                (float_of_int col_stats.max +. 1.) *. log (counts +. 1.)
              else begin
                if power < 1. then
                  (counts ** power) *. ((float_of_int col_stats.max) ** (1. -. power)) /. power
                else
                  (counts ** power)
              end in
            if quantize then
              floor v
            else
              v
    end
    (* Implementation function *)
    let stats_table_of_core_db ?(threads = 1) ?(verbose = false) core =
      let compute_one what n =
        let red_len =
          match what with
          | ColOrRow.Col ->
            core.n_rows - 1
          | Row ->
            core.n_cols - 1 in
        let non_zero = ref 0 and min = ref 0 and max = ref 0 and sum = ref 0. and sum_log = ref 0. in
        for i = 0 to red_len do
          let v =
            match what with
            | Col ->
              CountBAVector.N.to_int core.storage.(n).@(i)
            | Row ->
              CountBAVector.N.to_int core.storage.(i).@(n) in
          let f_v = float_of_int v in
          if f_v >= 0. then
            incr non_zero;
          min := Stdlib.min !min v;
          max := Stdlib.max !max v;
          sum := !sum +. f_v;
          sum_log := !sum_log +. log f_v
        done;
        { Transformation.non_zero = !non_zero;
          min = !min;
          max = !max;
          sum = !sum;
          sum_log = !sum_log } in
      let compute_all what =
        let n =
          match what with
          | ColOrRow.Col ->
            core.n_cols
          | Row ->
            core.n_rows in
        let step = n / threads / 5 |> max 1 and processed = ref 0 and res = ref [] in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if verbose then
              Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]%!"
                String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) !processed n;
            let to_do = n - !processed in
            if to_do > 0 then begin
              let to_do = min to_do step in
              let res = !processed, to_do in
              processed := !processed + to_do;
              res
            end else
              raise End_of_file)
          (fun (processed, to_do) ->
            let res = ref [] in
            for i = 0 to to_do - 1 do
              processed + i |> compute_one what |> List.accum res
            done;
            Array.of_rlist !res)
          (List.accum res)
          threads;
        let rec binary_merge_arrays processed to_do =
          match processed, to_do with
          | [], [] ->
            [||]
          | [], [a] ->
            a
          | [res], [] ->
            res
          | [res], [a] ->
            Array.append res a
          | _, [] ->
            List.rev processed |> binary_merge_arrays []
          | _, [a] ->
            a :: processed |> List.rev |> binary_merge_arrays []
          | _, a1 :: a2 :: tl ->
            binary_merge_arrays ((Array.append a1 a2) :: processed) tl in
        let res = List.rev !res |> binary_merge_arrays [] in
        if verbose then
          Printf.eprintf "%s\r(%s): Computing %s statistics [%d/%d]\n%!"
            String.TermIO.clear __FUNCTION__ (ColOrRow.to_string what) n n;
        res in
      { Transformation.col_stats = compute_all Col;
        row_stats = compute_all Row }
    type t = {
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
    module CombinationCriterion =
      struct
        type t = RescaledMean | RescaledMedian
        let of_string = function
          | "mean" -> RescaledMean
          | "median" -> RescaledMedian
          | s ->
            Exception.raise_unrecognized_initializer __FUNCTION__ "Combination criterion" s
        let to_string = function
          | RescaledMean -> "mean"
          | RescaledMedian -> "median"
      end
    (* It should be OK to have the same label on both LHS and RHS,
        as a temporary vector is used to generate the combination *)
    let add_combined_selected ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db new_label selection criterion =
      let stats = stats_table_of_core_db ~threads ~verbose db.core and db = ref db in
      if verbose then
        Printf.eprintf "(%s): Adding/replacing spectrum '%s': [%!" __FUNCTION__ new_label;
      (* We allocate the result *)
      add_empty_column_if_needed db new_label;
      let new_col_idx = StringHashtbl.find !db.col_names_to_idx new_label in
      let new_col = !db.core.storage.(new_col_idx) in
      (* Computing valid labels and maximum normalisation across columns *)
      let found_cols = ref [] and max_norm = ref 0. in
      StringSet.iter
        (fun label ->
          if verbose then
            Printf.eprintf " '%s'%!" label;
          (* Some labels might be invalid *)
          match StringHashtbl.find_opt !db.col_names_to_idx label with
          | Some col_idx ->
            List.accum found_cols col_idx;
            max_norm := max !max_norm stats.col_stats.(col_idx).sum
          | None ->
            if verbose then
              Printf.eprintf "(NOT FOUND)%!")
        selection;
      let found_cols = Array.of_list !found_cols in (* We don't really care about the order here *)
      let num_found_cols = Array.length found_cols and max_norm = !max_norm and norm = ref 0. in
      if verbose then
        Printf.eprintf " ] n_found=%d max_norm=%.16g.\n%!" num_found_cols max_norm;
      let n_rows = !db.core.n_rows and processed_rows = ref 0 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_rows < n_rows then
            let to_do = max 1 (elements_per_step / num_found_cols) |> min (n_rows - !processed_rows) in
            let new_processed_rows = !processed_rows + to_do in
            let res = !processed_rows, new_processed_rows - 1 in
            processed_rows := new_processed_rows;
            res
          else
            raise End_of_file)
        (fun (lo_row, hi_row) ->
          let module Freqs = Numbers.FloatFreqsVector in
          (* We need one histogram per row to be able to compute statistics such as the median *)
          let row_combinators = Array.init (hi_row - lo_row + 1) (fun _ -> Freqs.make ~non_negative:true ()) in
          for i = lo_row to hi_row do
            (* We iterate over valid columns *)
            Array.iter
              (fun col_idx ->
                (* We normalise columns *separately* before combining them *)
                let col = !db.core.storage.(col_idx) and norm = stats.col_stats.(col_idx).sum in
                (* All counts are non-negative *)
                if norm > 0. then
                  (* We add the renormalised sum to the suitable row histogram *)
                  CountBAVector.N.to_float col.@(i) *. max_norm /. norm |> Freqs.add row_combinators.(i - lo_row))
              found_cols
          done;
          (* For each row histogram in the input range, we now generate a combination and pass it along *)
          lo_row,
          Array.map
            (fun combinator ->
              match criterion with
              | CombinationCriterion.RescaledMean ->
                Freqs.sum combinator
              | RescaledMedian ->
                Freqs.median combinator *. float_of_int num_found_cols)
            row_combinators)
        (fun (lo_row, block) ->
          let n_processed = Array.length block in
          for i = lo_row to lo_row + n_processed - 1 do
            let res_i = block.(i - lo_row) in
            Numbers.Float.(norm ++ res_i);
            let res_i = CountBAVector.N.of_float res_i in
            (* Actual copy to storage *)
            new_col.@(i) <- res_i
          done;
          let old_processed_rows = !processed_rows in
          processed_rows := !processed_rows + n_processed;
          if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
            Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines%!"
              String.TermIO.clear __FUNCTION__ !processed_rows n_rows)
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Combining spectra: done %d/%d lines. Norm=%.16g\n%!"
          String.TermIO.clear __FUNCTION__ !processed_rows n_rows !norm;
      (* If metadata is present in the database, we generate some for the new column too *)
      if !db.core.n_meta > 0 then begin
        (* For each metadata field, we compute the intersection of the values across all selected columns *)
        let res = Array.make !db.core.n_meta StringSet.empty in
        StringSet.iter
          (fun label ->
            match StringHashtbl.find_opt !db.col_names_to_idx label with
            | Some col_idx ->
              let col = !db.core.meta.(col_idx) in
              for i = 0 to !db.core.n_meta - 1 do
                res.(i) <- StringSet.add col.(i) res.(i)
              done
            | None -> ())
          selection;
        !db.core.meta.(new_col_idx) <-
          Array.map
            (fun set ->
              if StringSet.cardinal set = 1 then
                StringSet.min_elt set
              else
                "")
            res
      end;
      !db
    exception Classes_label_not_found of string
    let get_indicator_vector db classes_label =
      match StringHashtbl.find_opt db.meta_names_to_idx classes_label with
      | None ->
        Classes_label_not_found classes_label |> raise
      | Some classes_label_idx ->
        (* We derive the class indicator vector *)
        let n_samples = db.core.n_cols in
        let num_classes = ref 0 and class_to_ind = ref StringMap.empty and ind_to_class = ref IntMap.empty
        and res = Array.make n_samples (-1) in
        (* We iterate explicitly as there might be trailing space at the end of vectors *)
        for i = 0 to n_samples - 1 do
          res.(i) <- begin
            let class_label = db.core.meta.(i).(classes_label_idx) in
            match StringMap.find_opt class_label !class_to_ind with
            | None ->
              (* We insert the new label *)
              let res = !num_classes in
              class_to_ind := StringMap.add class_label res !class_to_ind;
              ind_to_class := IntMap.add res class_label !ind_to_class;
              incr num_classes;
              res
            | Some ind ->
              ind
          end
        done;
        !num_classes, !ind_to_class, res
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
    exception Class_label_is_also_spectrum_name of string
    let split_spectra ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db classes_label criterion =
      let _, ind_to_class, ind_classes = get_indicator_vector db classes_label in
      if verbose then begin
        Printf.eprintf "(%s): Classes=[" __FUNCTION__;
        Array.iter (Printf.eprintf " %d") ind_classes;
        Printf.eprintf " ]\n%!"
      end;
      (* We split spectra names according to their class *)
      let module I2SMM = Tools.Multimap (ComparableInt) (ComparableString) in
      let split_names = ref I2SMM.empty in
      Array.iteri
        (fun i ind ->
          split_names := I2SMM.add ind db.core.idx_to_col_names.(i) !split_names)
        ind_classes;
      let res = ref db in
      I2SMM.iter_set
        (fun ind spectra_names ->
          let class_name = IntMap.find ind ind_to_class in
          if StringHashtbl.mem db.col_names_to_idx class_name then
            Class_label_is_also_spectrum_name class_name |> raise;
          res := add_combined_selected ~threads ~elements_per_step ~verbose !res class_name spectra_names criterion)
        !split_names;
      remove_selected !res (Array.to_seq db.core.idx_to_col_names |> StringSet.of_seq)
    (* *)
    module Stats = Numbers.OnlineStats(Numbers.Float) (*(CountBAVector.N)*)
    module Freqs = Numbers.FloatFreqsVector
    module FAVector = Numbers.FAVector
    module LF_FAVector = Numbers.LinearFit(FAVector)
    exception Invalid_number_of_classes of int
    let distill_kmers ?(threads = 1) ?(elements_per_step = 10000) ?(verbose = false)
        db classes_label summary_prefix =
      let stats_table = stats_table_of_core_db ~threads ~verbose db.core in
      let n_classes, _, ind_classes = get_indicator_vector db classes_label and n_samples = db.core.n_cols in
      if n_classes = 1 || n_classes = n_samples then
        Invalid_number_of_classes n_classes |> raise;
      if verbose then begin
        Printf.eprintf "(%s): Classes=[" __FUNCTION__;
        Array.iter (Printf.eprintf " %d") ind_classes;
        Printf.eprintf " ]\n%!"
      end;
      let n_kmers = db.core.n_rows and processed_kmers = ref 0 in
      let stats_classes = Array.init n_classes (fun _ -> Array.init n_classes (fun _ -> Stats.make ()))
      and avgs_on_mean = FAVector.(make n_kmers N.zero) and avgs_off_mean = FAVector.(make n_kmers N.zero)
      and avgs_on_median = FAVector.(make n_kmers N.zero) and avgs_off_median = FAVector.(make n_kmers N.zero)
      and vars_on_mean = FAVector.(make n_kmers N.zero) and vars_off_mean = FAVector.(make n_kmers N.zero)
      and vars_on_median = FAVector.(make n_kmers N.zero) and vars_off_median = FAVector.(make n_kmers N.zero)
      and covs_on_mean = FAVector.(make n_kmers N.zero) and covs_off_mean = FAVector.(make n_kmers N.zero)
      and covs_on_median = FAVector.(make n_kmers N.zero) and covs_off_median = FAVector.(make n_kmers N.zero) in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_kmers < n_kmers then
            let to_do =
              max 1 (elements_per_step / (n_samples * (n_samples - 1) / 2)) |> min (n_kmers - !processed_kmers) in
            let new_processed_kmers = !processed_kmers + to_do in
            let res = !processed_kmers, new_processed_kmers - 1 in
            processed_kmers := new_processed_kmers;
            res
          else
            raise End_of_file)
        (fun (lo_kmer, hi_kmer) ->
          let res = ref [] in
          for kmer = lo_kmer to hi_kmer do
            CountBAVector.(
              (* We iterate over all the couples of samples *)
              for i = 0 to n_samples - 1 do
                let ind_class_i = ind_classes.(i)
                (* Counts must be normalised *)
                and counts_i = N.to_float db.core.storage.(i).@(kmer) /. stats_table.col_stats.(i).sum in
                for j = i + 1 to n_samples - 1 do
                  let ind_class_j = ind_classes.(j)
                  (* Counts must be normalised *)
                  and counts_j = N.to_float db.core.storage.(j).@(kmer) /. stats_table.col_stats.(j).sum in
                  let diff = Float.abs (counts_i -. counts_j) in
                  if ind_class_i = ind_class_j then
                    Stats.add stats_classes.(ind_class_i).(ind_class_j) diff
                  else begin
                    if ind_class_i < ind_class_j then
                      Stats.add stats_classes.(ind_class_i).(ind_class_j) diff
                    else
                      Stats.add stats_classes.(ind_class_j).(ind_class_i) diff
                  end
                done
              done;
              let avgs_on = Freqs.make () and avgs_off = Freqs.make ()
              and vars_on = Freqs.make () and vars_off = Freqs.make ()
              and covs_on = Freqs.make () and covs_off = Freqs.make () in
              (* We iterate over all the couples of classes *)
              for i = 0 to n_classes - 1 do
                let stats = stats_classes.(i).(i) in
                Stats.mean stats |> Freqs.add avgs_on;
                Stats.sample_variance stats |> Freqs.add vars_on;
                Stats.sample_coefficient_of_variation stats |> Freqs.add covs_on;
                Stats.clear stats;
                for j = i + 1 to n_classes - 1 do
                  let stats = stats_classes.(i).(j) in
                  Stats.mean stats |> Freqs.add avgs_off;
                  Stats.sample_variance stats |> Freqs.add vars_off;
                  Stats.sample_coefficient_of_variation stats |> Freqs.add covs_off;
                  Stats.clear stats
                done
              done;
              begin
                kmer, Freqs.mean avgs_on, Freqs.mean avgs_off, Freqs.median avgs_on, Freqs.median avgs_off,
                      Freqs.mean vars_on, Freqs.mean vars_off, Freqs.median vars_on, Freqs.median vars_off,
                      Freqs.mean covs_on, Freqs.mean covs_off, Freqs.median covs_on, Freqs.median covs_off
              end |> List.accum res
            );
          done;
          !res)
        (fun block ->
          let old_processed_kmers = !processed_kmers in
          List.iter
            (fun
              (kmer,
               avgs_on_mean_kmer, avgs_off_mean_kmer, avgs_on_median_kmer, avgs_off_median_kmer,
               vars_on_mean_kmer, vars_off_mean_kmer, vars_on_median_kmer, vars_off_median_kmer,
               covs_on_mean_kmer, covs_off_mean_kmer, covs_on_median_kmer, covs_off_median_kmer) ->
              FAVector.(
                avgs_on_mean.@(kmer) <- avgs_on_mean_kmer;
                avgs_off_mean.@(kmer) <- avgs_off_mean_kmer;
                avgs_on_median.@(kmer) <- avgs_on_median_kmer;
                avgs_off_median.@(kmer) <- avgs_off_median_kmer;
                vars_on_mean.@(kmer) <- vars_on_mean_kmer;
                vars_off_mean.@(kmer) <- vars_off_mean_kmer;
                vars_on_median.@(kmer) <- vars_on_median_kmer;
                vars_off_median.@(kmer) <- vars_off_median_kmer;
                covs_on_mean.@(kmer) <- covs_on_mean_kmer;
                covs_off_mean.@(kmer) <- covs_off_mean_kmer;
                covs_on_median.@(kmer) <- covs_on_median_kmer;
                covs_off_median.@(kmer) <- covs_off_median_kmer
              );
              incr processed_kmers
          ) block;
          if verbose && !processed_kmers / 100 > old_processed_kmers / 100 then
            Printf.eprintf "%s\r(%s): Distilled %d/%d kmers%!"
              String.TermIO.clear __FUNCTION__ !processed_kmers n_kmers)
        threads;
      let fit_avgs_mean, _, residuals_avgs_mean = LF_FAVector.make avgs_on_mean avgs_off_mean
      and fit_avgs_median, _, residuals_avgs_median = LF_FAVector.make avgs_on_median avgs_off_median
      and fit_vars_mean, _, residuals_vars_mean = LF_FAVector.make vars_on_mean vars_off_mean
      and fit_vars_median, _, residuals_vars_median = LF_FAVector.make vars_on_median vars_off_median
      and fit_covs_mean, _, residuals_covs_mean = LF_FAVector.make covs_on_mean covs_off_mean
      and fit_covs_median, _, residuals_covs_median = LF_FAVector.make covs_on_median covs_off_median in
      if verbose then begin
        Printf.eprintf "%s\r(%s): Distilled %d/%d kmers.\n"
          String.TermIO.clear __FUNCTION__ !processed_kmers n_kmers;
        Printf.eprintf "(%s): Fit for avgs mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_avgs_mean) (LF_FAVector.get_slope fit_avgs_mean);
        Printf.eprintf "(%s): Fit for avgs median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_avgs_median) (LF_FAVector.get_slope fit_avgs_median);
        Printf.eprintf "(%s): Fit for vars mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_vars_mean) (LF_FAVector.get_slope fit_vars_mean);
        Printf.eprintf "(%s): Fit for vars median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_vars_median) (LF_FAVector.get_slope fit_vars_median);
        Printf.eprintf "(%s): Fit for covs mean is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_covs_mean) (LF_FAVector.get_slope fit_covs_mean);
        Printf.eprintf "(%s): Fit for covs median is %.6g + %.6g * x\n%!" __FUNCTION__
          (LF_FAVector.get_intercept fit_covs_median) (LF_FAVector.get_slope fit_covs_median)
      end;
      (* We output the summary *)
      let summary = {
        Matrix.which = Distill;
        matrix = {
          col_names = Array.sub db.core.idx_to_row_names 0 db.core.n_rows;
          row_names =
            [| "InnerAvgMean"; "OuterAvgMean"; "ResidualAvgMean";
               "InnerAvgMedian"; "OuterAvgMedian"; "ResidualAvgMedian";
               "InnerVarMean"; "OuterVarMean"; "ResidualVarMean";
               "InnerVarMedian"; "OuterVarMedian"; "ResidualVarMedian";
               "InnerCOVMean"; "OuterCOVMean"; "ResidualCOVMean";
               "InnerCOVMedian"; "OuterCOVMedian"; "ResidualCOVMedian" |];
          data = FAVector.(
            [| to_floatarray avgs_on_mean; to_floatarray avgs_off_mean; to_floatarray residuals_avgs_mean;
               to_floatarray avgs_on_median; to_floatarray avgs_off_median; to_floatarray residuals_avgs_median;
               to_floatarray vars_on_mean; to_floatarray vars_off_mean; to_floatarray residuals_vars_mean;
               to_floatarray vars_on_median; to_floatarray vars_off_median; to_floatarray residuals_vars_median;
               to_floatarray covs_on_mean; to_floatarray covs_off_mean; to_floatarray residuals_covs_mean;
               to_floatarray covs_on_median; to_floatarray covs_off_median; to_floatarray residuals_covs_median |]
          )
        }
      } in
      Matrix.to_file ~threads (Matrix.transpose ~threads summary) summary_prefix
    let transform ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) transformation db =
      let transform = Transformation.compute ~which:transformation
      and stats = stats_table_of_core_db ~threads ~verbose db.core
      and n_rows = db.core.n_rows and n_cols = db.core.n_cols
      and new_storage = Array.copy db.core.storage and processed_cols = ref 0 in
      Processes.Parallel.process_stream_chunkwise
        (fun () ->
          if !processed_cols < n_cols then
            let to_do = max 1 (elements_per_step / n_rows) |> min (n_cols - !processed_cols) in
            let new_processed_cols = !processed_cols + to_do in
            let res = !processed_cols, new_processed_cols - 1 in
            processed_cols := new_processed_cols;
            res
          else
            raise End_of_file)
        (fun (lo_col, hi_col) ->
          let res = ref [] in
          for col_idx = lo_col to hi_col do
            List.accum res begin
              col_idx,
              CountBAVector.mapi
                (fun row_idx n ->
                  transform
                    ~col_num:n_cols ~col_stats:stats.col_stats.(col_idx)
                    ~row_num:n_rows ~row_stats:stats.row_stats.(row_idx) n)
                db.core.storage.(col_idx)
            end
          done;
          !res)
        (List.iter
          (fun (col_idx, col) ->
            new_storage.(col_idx) <- col;
            incr processed_cols;
            if verbose then
              Printf.eprintf "%s\r(%s): Transforming database: done %d/%d %s%!"
                String.TermIO.clear __FUNCTION__ !processed_cols n_cols
                (String.pluralize_int "column" !processed_cols)))
        threads;
      if verbose then
        Printf.eprintf "%s\r(%s): Transforming database: done %d/%d %s.\n%!"
          String.TermIO.clear __FUNCTION__ n_cols n_cols
          (String.pluralize_int "column" n_cols);
      { db with core = { db.core with storage = new_storage } }
    let to_distances
        ?(precision = 15) ?(normalise = true) ?(threads = 1) ?(elements_per_step = 100) ?(verbose = false)
        distance db selection_1 selection_2 prefix =
      let stats = stats_table_of_core_db ~threads ~verbose db.core
      and n_r = db.core.n_rows in (* Number of k-mers. It stays the same even if we select a subset of spectra *)
      let make_submatrix selection =
        (* We select spectra, i.e. columns *)
        let idxs = ref IntSet.empty in
        Array.iteri
          (fun col_idx col_name ->
            if StringSet.mem col_name selection then
              idxs := IntSet.add col_idx !idxs)
          db.core.idx_to_col_names;
        let idxs = IntSet.elements_array !idxs in
        let n_c = Array.length idxs in (* We now know the number of columns as well *)
        (* The distance matrix is computed rowwise, and k-mers are physically stored as rows in db:
            we need to transpose *)
        { Matrix.Base.col_names = Array.sub db.core.idx_to_row_names 0 n_r;
          row_names = Array.init n_c (fun i -> db.core.idx_to_col_names.(idxs.(i)));
          data =
            (* Here we just need to convert the counts to floats, and possibly normalise *)
            Array.init n_c
              (fun i ->
                let idx = idxs.(i) in
                let norm = stats.col_stats.(idx).sum in
                let norm =
                  if not normalise || norm = 0. then
                    1.
                  else
                    norm in
                (*
                Printf.eprintf "Norm@%d=%g\n%!" idx norm;
                *)
                Float.Array.init n_r
                  (fun j -> CountBAVector.N.to_float db.core.storage.(idx).@(j) /. norm)) } in
      let metric = Float.Array.make n_r 1.
      and matrix_1 = make_submatrix selection_1 and matrix_2 = make_submatrix selection_2 in
      Matrix.to_file ~precision ~threads ~elements_per_step ~verbose {
        which = DMatrix;
        matrix =
          Matrix.Base.get_distance_rowwise ~threads ~elements_per_step ~verbose distance metric matrix_1 matrix_2
      } prefix
    (* *)
    let make_filename_table = function
      | w when String.length w >= 5 && String.sub w 0 5 = "/dev/" -> w
      | prefix -> prefix ^ ".KPopCounter.txt"
    let to_files ?(precision = 15) ?(output_zero_kmers = false)
                 ?(threads = 1) ?(elements_per_step = 40000) ?(verbose = false) db prefix =
      let stats = stats_table_of_core_db ~threads ~verbose db.core
      and fname = make_filename_table prefix in
      let output = open_out fname and meta = ref [] and rows = ref [] and cols = ref [] in
      (* We determine which rows and colunms should be output after all filters have been applied *)
      (*  Rows: metadata and k-mers *)
      Array.iteri
        (fun i meta_name ->
          (* There might be additional storage *)
          if i < db.core.n_meta then
            List.accum meta (meta_name, i))
        db.core.idx_to_meta_names;
      Array.iteri
        (fun i row_name ->
          (* There might be additional storage.
             Also, we only print the row if it has non-zero elements or if we are explicitly requested to do so *)
          if i < db.core.n_rows && (stats.row_stats.(i).sum > 0. || output_zero_kmers) then
            List.accum rows (row_name, i))
        db.core.idx_to_row_names;
      (*  Columns *)
      Array.iteri
        (fun i col_name ->
          (* There might be additional storage *)
          if i < db.core.n_cols then
            List.accum cols (col_name, i))
        db.core.idx_to_col_names;
      let meta = Array.of_rlist !meta and rows = Array.of_rlist !rows and cols = Array.of_rlist !cols in
      let n_rows = Array.length rows and n_cols = Array.length cols in
      (* There must be at least one row to print.
         If there are no columns, we just print metadata/row names  *)
      if (Array.length meta + n_rows) > 0 then begin
        (* We print column names *)
        Array.iter
          (fun (col_name, _) -> Printf.fprintf output "\t%s" col_name)
          cols;
        Printf.fprintf output "\n";
        (* We print metadata as lines *)
        Array.iter
          (fun (meta_name, meta_idx) ->
            Printf.fprintf output "%s" meta_name;
            Array.iter
              (fun (_, col_idx) -> Printf.fprintf output "\t%s" db.core.meta.(col_idx).(meta_idx))
              cols;
            Printf.fprintf output "\n")
          meta;
        Printf.fprintf output "%!";
        let processed_rows = ref 0 and buf = Buffer.create 1048576 in
        Processes.Parallel.process_stream_chunkwise
          (fun () ->
            if !processed_rows < n_rows then
              let to_do = max 1 (elements_per_step / (max 1 n_cols)) |> min (n_rows - !processed_rows) in
              let new_processed_rows = !processed_rows + to_do in
              let res = !processed_rows, new_processed_rows - 1 in
              processed_rows := new_processed_rows;
              res
            else
              raise End_of_file)
          (fun (lo_row, hi_row) ->
            Buffer.clear buf;
            for i = lo_row to hi_row do
              let row_name, row_idx = rows.(i) in
              Printf.bprintf buf "%s" row_name;
              Array.iter
                (fun (_, col_idx) ->
                  Printf.bprintf buf "\t%.*g" precision db.core.storage.(col_idx).@(row_idx))
                cols;
              Printf.bprintf buf "\n"
            done;
            hi_row - lo_row + 1, Buffer.contents buf)
          (fun (n_processed, block) ->
            Printf.fprintf output "%s" block;
            let old_processed_rows = !processed_rows in
            processed_rows := !processed_rows + n_processed;
            if verbose && !processed_rows / 10000 > old_processed_rows / 10000 then
              Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d lines%!"
                String.TermIO.clear __FUNCTION__ fname !processed_rows n_rows)
          threads;
        if verbose then
          Printf.eprintf "%s\r(%s): Writing table to file '%s': done %d/%d lines.\n%!"
            String.TermIO.clear __FUNCTION__ fname n_rows n_rows
      end;
      close_out output
    let of_files ?(threads = 1) ?(bytes_per_step = 4194304) ?(verbose = false) prefix =


      make_empty ()


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
    module Transformation:
      sig
        type statistics_t = {
          non_zero: int;
          min: int;
          max: int;
          sum: float;
          sum_log: float;
          (*
          mean: float;
          variance: float
          (* Could be extended with more moments if needed *)
          harmonic: float (* The harmonic mean of non-zero elements *)
          *)
        }
        (* Transformation function *)
        type t =
          (* Thresholds in the interval (0.,1.) are taken as relative ones *)
          | Threshold of float
          | Power of float
          | CLR
          | Pseudo of float * bool (* The boolean is whether to quantise *)
        val of_string: string -> t
        val to_string: t -> string
      end
    type t = {
      core: marshalled_t;
      (* Inverted hashes for parsing *)
      col_names_to_idx: int StringHashtbl.t; (* Labels *)
      row_names_to_idx: int StringHashtbl.t; (* Hashes *)
      meta_names_to_idx: int StringHashtbl.t (* Metadata fields *)
    }
    val make_empty: unit -> t
    val merge: ?verbose:bool -> t -> t -> t
    (* Increase the counts for the specified sample and k-mer by the given amount *)
    val add_counts: t -> string -> string -> int -> t
    (* Add metadata - the first field must be the label *)
    exception Wrong_number_of_columns of int * int * int
    val add_metadata: ?verbose:bool -> t -> string -> t
    (* Select column names identified by regexps on metadata fields *)
    val selected_from_regexps: ?verbose:bool -> t -> (string * Str.regexp) list -> StringSet.t
    val selected_negate: t -> StringSet.t -> StringSet.t
    module CombinationCriterion:
      sig
        type t = RescaledMean | RescaledMedian
        val of_string: string -> t
        val to_string: t -> string
      end
    (* Generate according to the specified criterion a combination of the spectra having the given labels,
        name the combination as directed, and add it to the database *)
    val add_combined_selected: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                               t -> string -> StringSet.t -> CombinationCriterion.t -> t
    (* Combine spectra into class representatives according to the specified class labels *)
    exception Class_label_is_also_spectrum_name of string
    val split_spectra: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                       t -> string -> CombinationCriterion.t -> t
    (* Remove spectra with the given labels *)
    val remove_selected: t -> StringSet.t -> t
    (* Distill k-mers, i.e., sort them by decreasing discriminative power according to the specified class labels *)
    exception Classes_label_not_found of string
    exception Invalid_number_of_classes of int
    val distill_kmers: ?threads:int -> ?elements_per_step:int -> ?verbose:bool -> t -> string -> string -> unit
    (* Transformations *)
    val transform: ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                   Transformation.t -> t -> t
    (* Spectral distance matrix *)
    val to_distances: ?precision:int -> ?normalise:bool -> ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                      Space.Distance.t -> t -> StringSet.t -> StringSet.t -> string -> unit
    (* Output information about the contents *)
    val output_summary: ?verbose:bool -> t -> unit
    (* *)
    val to_files: ?precision:int -> ?output_zero_kmers:bool ->
                  ?threads:int -> ?elements_per_step:int -> ?verbose:bool ->
                  t -> string -> unit
    val of_files: ?threads:int -> ?bytes_per_step:int -> ?verbose:bool -> string -> t
    val to_binary: ?verbose:bool -> t -> string -> unit
    val of_binary: ?verbose:bool -> string -> t (* Can fail due to archive version *)
  end
)

