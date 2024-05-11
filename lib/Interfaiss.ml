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

include (
  struct
    type index_t =
      | Flat
      | PQ of int * int
      | HNSW of int * int
    type t
    external create: index_t -> int -> t = "InterfaissCreate"
    let create ?(index_type = Flat) n =
      try
        create index_type n
      with _ ->
        (* This might happen if the definition of index_t get out of synchro with the C interface *)
        assert false
    type vectors_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    type offsets_t = (int, Bigarray.int64_elt, Bigarray.c_layout) Bigarray.Array2.t
    external add: t -> vectors_t -> unit = "InterfaissAdd"
    external train: t -> vectors_t -> unit = "InterfaissTrain"
    external query: t -> vectors_t -> int -> vectors_t * offsets_t = "InterfaissQuery"
    external delete: t -> unit = "InterfaissDelete"
  end: sig
    type index_t =
      | Flat
      | PQ of int * int
      | HNSW of int * int
    type t
    val create: ?index_type:index_t -> int -> t
    type vectors_t = (float, Bigarray.float32_elt, Bigarray.c_layout) Bigarray.Array2.t
    type offsets_t = (int, Bigarray.int64_elt, Bigarray.c_layout) Bigarray.Array2.t
    val add: t -> vectors_t -> unit
    val train: t -> vectors_t -> unit
    val query: t -> vectors_t -> int -> vectors_t * offsets_t
    val delete: t -> unit
  end
)

