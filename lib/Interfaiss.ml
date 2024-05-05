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

include (
  struct
    module F32BAVector = Numbers.Bigarray.Vector (
      struct
        include Numbers.Float
        type elt_t = Bigarray.float32_elt
        let elt = Bigarray.Float32
      end
    )
    type index_t =
      | Flat
      | PQ of int * int
      | HSNW of int * int
    type t
    external make: index_t -> F32BAVector.t -> t = "InterfaissMake"
    let make ?(index_type = Flat) vectors = make index_type vectors

    module XYZ = IntMap

  end (*: sig
    module F32BAVector = Numbers.Bigarray.Vector (
      struct
        include Numbers.Float
        type elt_t = Bigarray.float32_elt
        let elt = Bigarray.Float32
      end
    )
    type t
    type index_t =
      | Flat
      | PQ of int * int
      | HSNW of int * int
    val make ?index_type:index_t -> Bigarray -> t
  end *)
)


